#!/bin/bash
# 断点续传: 对 .parts_* 目录里未完成的 p<seg> 段, 按 z<seg>_<j> 子文件逐个补齐尾部
# 用法: _resume.sh <URL> <总字节> <输出文件> <线程数> <分段数N> <起始段号> <结束段号>
set -u
URL="$1"; TOTAL="$2"; OUT="$3"; SUB="${4:-6}"; N="${5:-28}"; A="${6:-0}"; B="${7:-27}"

WORK="$(dirname "$OUT")/.parts_$(basename "$OUT")"
mkdir -p "$WORK"
SEG=$(( (TOTAL + N - 1) / N ))

for ((seg=A; seg<=B; seg++)); do
  s=$(( seg*SEG )); [ "$s" -ge "$TOTAL" ] && continue
  e=$(( s + SEG - 1 )); [ "$e" -ge "$TOTAL" ] && e=$((TOTAL-1))
  SS=$(( (e - s + 1 + SUB - 1) / SUB ))
  echo "--- p$seg: $s-$e (want $((e-s+1)), SS=$SS)"
  ok=1
  for ((j=0; j<SUB; j++)); do
    ss=$(( s + j*SS )); [ "$ss" -gt "$e" ] && continue
    ee=$(( ss + SS - 1 )); [ "$ee" -gt "$e" ] && ee=$e
    want=$(( ee - ss + 1 ))
    f="$WORK/z${seg}_$j"
    [ -f "$WORK/p$seg" ] && { echo "  p$seg 已存在, 跳过"; ok=1; break; }
    have=0; [ -f "$f" ] && have=$(stat -c %s "$f")
    if [ "$have" -ge "$want" ]; then
      [ "$have" -ne "$want" ] && { echo "  z${seg}_$j 超长 $have>$want, 截断"; truncate -s "$want" "$f"; }
      continue
    fi
    # 只补尾部
    ns=$(( ss + have ))
    curl -s --retry 6 --retry-delay 3 --retry-all-errors -m 1800 \
         -r "$ns-$ee" -o "$WORK/.tmp_${seg}_$j" "$URL" &
  done
  wait
  for ((j=0; j<SUB; j++)); do
    ss=$(( s + j*SS )); [ "$ss" -gt "$e" ] && continue
    ee=$(( ss + SS - 1 )); [ "$ee" -gt "$e" ] && ee=$e
    want=$(( ee - ss + 1 ))
    f="$WORK/z${seg}_$j"
    [ -f "$WORK/.tmp_${seg}_$j" ] && { cat "$WORK/.tmp_${seg}_$j" >> "$f"; rm -f "$WORK/.tmp_${seg}_$j"; }
    have=0; [ -f "$f" ] && have=$(stat -c %s "$f")
    if [ "$have" -ne "$want" ]; then echo "  z${seg}_$j 仍缺: $have/$want"; ok=0; else echo "  z${seg}_$j OK"; fi
  done
  if [ "$ok" -eq 1 ]; then
    : > "$WORK/p$seg"
    for ((j=0; j<SUB; j++)); do [ -f "$WORK/z${seg}_$j" ] && cat "$WORK/z${seg}_$j" >> "$WORK/p$seg"; done
    got=$(stat -c %s "$WORK/p$seg")
    echo "  p$seg 合成 $got / $((e-s+1))"
    [ "$got" -ne "$((e-s+1))" ] && { echo "  p$seg 大小不符!"; rm -f "$WORK/p$seg"; }
  fi
done
echo "续传结束"
