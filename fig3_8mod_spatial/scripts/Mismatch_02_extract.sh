#!/bin/bash
# A->G mismatch extraction driver: stream modkit bed.gz, capture target sites
# + background A-position sample (Nvalid>=15, p=0.003).
# 2026-10-22. ~2min per 36M-row file over USB.
set -u
BASE=/mnt/f/BN2505012AH02S66N1/DirectRNA/05.methy
TG=/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/rediportal_20261022/sitelevel/mismatch/target_sites.bedkeys.tsv
OUTD=/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/rediportal_20261022/sitelevel/mismatch
mkdir -p "$OUTD"

extract() {  # $1=mod dir, $2=sample
  local bed="$BASE/$1/methy_detect/$2/$2.$1.bed.gz"
  local out="$OUTD/extract_$1_$2.tsv"
  [ -f "$bed" ] || { echo "MISSING $bed"; return 1; }
  zcat "$bed" | awk -F'\t' -v tg="$TG" -v bgp=0.003 -v cls_filter="$1" '
    BEGIN{while((getline line < tg)>0){split(line,a,"\t"); T[a[1] SUBSEP a[2]]=a[3]}
          srand(42)}
    {k=$1 SUBSEP $2
     if(k in T){ c=T[k]; if(c ~ /^ino_/ && cls_filter=="inosine") print c"\t"$0;
                 else if(c=="m6A_ctrl" && cls_filter=="m6a") print c"\t"$0 }
     else if($10>=15 && rand()<bgp) print "bg\t"$0 }' > "$out"
  echo "DONE $1/$2 -> $(wc -l < "$out") rows"
}

for s in NCs1 NCs2 NCs3 Day1s1 Day1s2 Day1s3; do
  extract inosine "$s"
done
for s in NCs1 NCs2 NCs3 Day1s1 Day1s2 Day1s3; do
  extract m6a "$s"
done
echo "ALL EXTRACTIONS DONE"
