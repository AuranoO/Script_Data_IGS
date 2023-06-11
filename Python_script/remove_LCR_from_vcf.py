with open("../ref_kamtchatka/Cedratvirus_kamchatka_longread_masked_position", "r") as maskpos:
    pos = maskpos.readlines()
list_interval = []
for line in pos:
    if line.startswith('>'):
        pass
    else:
        line = line.replace("\n", "")
        line = line.split(" ")
        interval = [line[0], line[2]]
        list_interval.append(interval)

print(list_interval)

with open("../variant_calling/lofreq_rnaseq2_filter.vcf", "r") as vcf:
    vcf_read = vcf.readlines()

snp_to_keep = open("../variant_calling/lofreq_rnaseq2_filter_noLCR.vcf", "w")
count = 0
for line in vcf_read:
    if not line.startswith("#"):
        split = line.split("\t")
        position = split[1]
        keep = 0
        for interval in list_interval:
            left = interval[0]
            right = interval[1]
            if left <= position <= right:
                count += 1
                keep = 1
        if keep == 0:
            snp_to_keep.write(line)
            pass
    else:
        snp_to_keep.write(line)
        pass
print(count)
snp_to_keep.close()

