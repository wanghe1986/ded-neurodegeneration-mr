# -*- coding: utf-8 -*-
"""eQTLGen 数据并发分段下载（单连接实测仅 ~9 KB/s）"""
import urllib.request, ssl, os, time, sys, threading
from concurrent.futures import ThreadPoolExecutor

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0 (research)"}
BASE = "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"

JOBS = [
    ("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
     "data/raw/eqtlgen_cis_fdr005.txt.gz"),
    ("2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz",
     "data/raw/eqtlgen_maf.txt.gz"),
]
NTHREAD = 16


def head_size(url):
    r = urllib.request.urlopen(urllib.request.Request(url, headers=UA, method="HEAD"),
                               timeout=60, context=CTX)
    return int(r.headers["Content-Length"])


def fetch_range(url, start, end, dst, idx, retries=5):
    part = f"{dst}.part{idx:03d}"
    if os.path.exists(part):
        have = os.path.getsize(part)
        want = end - start + 1
        if have >= want:
            return want, 0
        start = start + have
        mode = "ab"
    else:
        mode = "wb"
    for a in range(retries):
        try:
            h = dict(UA); h["Range"] = f"bytes={start}-{end}"
            r = urllib.request.urlopen(urllib.request.Request(url, headers=h),
                                       timeout=180, context=CTX)
            got = 0
            with open(part, mode) as o:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    o.write(b); got += len(b)
            return os.path.getsize(part), got
        except Exception as e:
            if a == retries - 1:
                print(f"  [part{idx:03d}] 失败: {e}", flush=True)
                return os.path.getsize(part) if os.path.exists(part) else 0, 0
            time.sleep(2 + a * 3)
    return 0, 0


def download(url, dst, nthread=NTHREAD):
    total = head_size(url)
    print(f"{os.path.basename(dst)}: 总大小 {total/1e6:.1f} MB, {nthread} 线程", flush=True)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    chunk = total // nthread
    spans = []
    for i in range(nthread):
        s = i * chunk
        e = total - 1 if i == nthread - 1 else (i + 1) * chunk - 1
        spans.append((s, e, i))
    t0 = time.time()
    done = [0]
    lock = threading.Lock()

    def work(sp):
        s, e, i = sp
        _, got = fetch_range(url, s, e, dst, i)
        with lock:
            done[0] += 1
        return got

    with ThreadPoolExecutor(max_workers=nthread) as ex:
        list(ex.map(work, spans))

    # 合并分片
    with open(dst, "wb") as out:
        for i in range(nthread):
            p = f"{dst}.part{i:03d}"
            if os.path.exists(p):
                with open(p, "rb") as f:
                    while True:
                        b = f.read(1 << 22)
                        if not b:
                            break
                        out.write(b)
                os.remove(p)
    sz = os.path.getsize(dst)
    dt = time.time() - t0
    print(f"  完成 {dst}  {sz/1e6:.1f}/{total/1e6:.1f} MB  用时 {dt:.0f}s  "
          f"平均 {sz/1e6/max(dt,1):.2f} MB/s", flush=True)
    return sz == total


for fname, dst in JOBS:
    url = BASE + fname
    if os.path.exists(dst):
        print(f"已存在, 跳过: {dst} ({os.path.getsize(dst)/1e6:.1f} MB)", flush=True)
        continue
    ok = download(url, dst)
    print(f"  {'校验通过' if ok else '大小不符, 需重跑'}", flush=True)

print("ALL DONE", flush=True)
