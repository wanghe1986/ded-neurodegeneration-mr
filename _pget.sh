#!/bin/bash
# 多线程分段下载 (EBI FTP 支持 Range, 且是每连接限速 -> 并发越多越快)
# 用法: _pget.sh <URL> <总字节> <输出文件> <线程数> [已有前缀文件]
set -u
URL="$1"; TOTAL="$2"; OUT="$3"; N="${4:-24}"; PREFIX="${5:-}"

WORK="$(dirname "$OUT")/.parts_$(basename "$OUT")"
mkdir -p "$WORK"

SKIP=0
if [ -n "$PREFIX" ] && [ -f "$PREFIX" ]; then
  SKIP=$(stat -c %s "$PREFIX")
  cp "$PREFIX" "$WORK/prefix"
fi
echo "总大小 $TOTAL 字节, 已有前缀 $SKIP 字节, 待下载 $((TOTAL-SKIP)) 字节, 线程 $N"

REM=$((TOTAL - SKIP))
if [ "$REM" -le 0 ]; then
  cp "$WORK/prefix" "$OUT"; echo "已完整, 无需下载"; exit 0
fi
SEG=$(( (REM + N - 1) / N ))

pids=()
for ((i=0; i<N; i++)); do
  s=$((SKIP + i*SEG))
  [ "$s" -ge "$TOTAL" ] && continue
  e=$((s + SEG - 1)); [ "$e" -ge "$TOTAL" ] && e=$((TOTAL-1))
  curl -s --retry 5 --retry-delay 3 --retry-all-errors -m 3600 \
       -r "$s-$e" -o "$WORK/p$i" "$URL" &
  pids+=($!)
done
echo "已启动 ${#pids[@]} 个分段下载"
wait

# 校验每段
bad=0
for ((i=0; i<N; i++)); do
  s=$((SKIP + i*SEG)); [ "$s" -ge "$TOTAL" ] && continue
  e=$((s + SEG - 1)); [ "$e" -ge "$TOTAL" ] && e=$((TOTAL-1))
  want=$((e - s + 1))
  if [ ! -f "$WORK/p$i" ]; then echo "缺少段 p$i"; bad=1; continue; fi
  got=$(stat -c %s "$WORK/p$i")
  if [ "$got" -ne "$want" ]; then echo "段 p$i 不完整: $got / $want"; bad=1; fi
done
if [ "$bad" -ne 0 ]; then echo "存在不完整分段, 请重跑本脚本 (已完成的分段会被复用)"; exit 1; fi

# 合并
: > "$OUT"
[ -f "$WORK/prefix" ] && cat "$WORK/prefix" >> "$OUT"
for ((i=0; i<N; i++)); do [ -f "$WORK/p$i" ] && cat "$WORK/p$i" >> "$OUT"; done

got=$(stat -c %s "$OUT")
echo "合并完成: $got / $TOTAL"
if [ "$got" -eq "$TOTAL" ]; then
  echo "OK 大小校验通过"; rm -rf "$WORK"
else
  echo "大小不符!"; exit 1
fi
