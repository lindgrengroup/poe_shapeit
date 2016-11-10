source 0_config.sh

for i in {1..22}
do
  plink --file ${INPUT} --chr $i --maf 0.05 --geno 0.05 --out ${INPUT}.chr${i} --recode --list-duplicate-vars ids-only suppress-first
  plink --file ${INPUT}.chr${i} --exclude ${INPUT}.chr${i}.dupvar --out ${INPUT}.chr${i} --recode
  cut -f4 ${INPUT}.chr${i}.map | sort | uniq -d > ${INPUT}.chr${i}.exclude
done
