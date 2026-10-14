#!/bin/bash
# One-line status for every attack run still in flight. Counts come from DISK, never from a
# loop's own bookkeeping -- a filter cannot tell "found nothing" from "ran nothing".
L=/data/ds85/RanPanMixer/runs/beagle_attack/logs
R=/data/ds85/RanPanMixer/runs/beagle_attack/root
printf "%-26s %-12s %s\n" "run" "progress" "state"
printf "%-26s %-12s %s\n" "--------------------------" "------------" "-----"
g=$(ls $L/gap_score_row*.status 2>/dev/null | wc -l)
printf "%-26s %-12s %s\n" "gap_score" "$g/308" "DONE"
b=$(find $R -name accuracy_results_protocol.json 2>/dev/null | wc -l)
printf "%-26s %-12s %s\n" "Beagle (their code)" "$b/308" "DONE"
d=$(ls $L/diploid_gap_score_row*.status 2>/dev/null | wc -l)
printf "%-26s %-12s %s\n" "diploid_gap_score" "$d/308" "$(cat $L/DIP.status 2>/dev/null || echo running)"
m=$(ls $L/MIA_privacy_row*.status 2>/dev/null | wc -l)
printf "%-26s %-12s %s\n" "MIA_privacy" "$m/308" "$(cat $L/MIA.status 2>/dev/null || echo running)"
s=$(find $R -name accuracy_results_selfout.json 2>/dev/null | wc -l)
printf "%-26s %-12s %s\n" "Beagle (their paper)" "$s/308" "$(cat $L/SELFOUT.status 2>/dev/null || echo running)"
echo
f=$(cat $L/*_row*.status 2>/dev/null | grep -vc '^0$')
echo "failed jobs: $f"
[ "$f" -gt 0 ] && { echo "FAILURES:"; for x in $L/*_row*.status; do [ "$(cat $x)" != "0" ] && echo "  $(basename $x .status) -> $(cat $x)"; done; }
echo "live procs: diploid=$(pgrep -cf diploid_gap_score) mia=$(pgrep -cf MIA_privacy) beagle=$(pgrep -cf beagle_attack_one)"
uptime | sed 's/^/load: /'
