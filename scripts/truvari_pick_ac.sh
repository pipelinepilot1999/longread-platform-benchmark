# Sensitivity analysis: the 12 Truvari benchmarks again with --pick ac (a 1/1 call may match two truth records),
# as the GIAB v5.0q README recommends. Run on the host from anywhere; writes results/*/truvari_pick_ac/ and
# results/sv_benchmark_pick_ac.tsv. See docs/decisions.md D-024.
set -euo pipefail
REPO=/home/ec2-user/longread/repo; RES=/home/ec2-user/longread/resources
TR="docker run --rm --cpus 1 -u $(id -u):$(id -g) -v /home/ec2-user:/home/ec2-user -w $REPO longread-platform-benchmark/truvari:5.4.0 truvari"
cd $REPO
declare -A REG=( [permissive]="--refdist 1000 --pctseq 0 --pctsize 0.5 --typeignore" [default]="--refdist 500 --pctseq 0.7 --pctsize 0.7" [strict]="--refdist 100 --pctseq 0.9 --pctsize 0.9" )
printf "platform\tcaller\tregime\ttp_base\ttp_comp\tfp\tfn\tprecision\trecall\tf1\n" > results/sv_benchmark_pick_ac.tsv
for p in ont hifi; do for c in cutesv sniffles; do for r in permissive default strict; do
  out=results/$p/truvari_pick_ac/${c}_${r}; rm -rf $out; mkdir -p results/$p/truvari_pick_ac
  $TR bench -b $RES/truth_chr20/truth_sv.chr20.vcf.gz -c results/$p/sv/$p.$c.vcf.gz -f /home/ec2-user/wgs-wes/ref/GRCh38.fa \
     --includebed $RES/truth_chr20/truth_sv.chr20.bed --passonly --sizemin 50 --sizefilt 30 ${REG[$r]} --pick ac -o $out > /dev/null 2>&1
  echo "truvari bench --includebed truth_sv.chr20.bed --passonly --sizemin 50 --sizefilt 30 ${REG[$r]} --pick ac" > $out/cmd.txt
  python3 - "$p" "$c" "$r" "$out/summary.json" >> results/sv_benchmark_pick_ac.tsv <<'PY'
import json,sys
p,c,r,f=sys.argv[1:]; s=json.load(open(f))
print(p,c,r,s['TP-base'],s['TP-comp'],s['FP'],s['FN'],round(s['precision'],4),round(s['recall'],4),round(s['f1'],4),sep='\t')
PY
done; done; done
