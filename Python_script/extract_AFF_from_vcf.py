with open("../vcf2/lofreq_unicycler", "r") as vcf:
    vcf_read = vcf.readlines()

af_to_append = open("/home/toumi/Cedratvirus/vcf2/af_frequency.txt", "w")

# For LoFreq
for line in vcf_read:
    if line.startswith("#"):
        pass
    else:
        line = line.split("AF=")[1].split(';')
        af_to_append.write(line[0]+"\n")

# For FreeBayes
# a = 0
# for line in vcf_read:
#     if line.startswith("#"):
#         pass
#     else:
#         line = line.split("AB=")[1].split(';')[0]
#         line = line.split(",")
#         if len(line) > 1:
#             a += 1
#             line = max(line)
#             af_to_append.write(line + "\n")
#         else:
#             af_to_append.write(max(line) + "\n")
#
# print(a)
# af_to_append.close()

# variant_list = []
# for line in vcf_read:
#     if line.startswith("#"):
#         pass
#     else:
#         line = line.split("EFF=")[1].split("(")[1].split("|")
#         print(line)
#         variant_list.append(line[1])
#
# unique = set(variant_list)
# print(unique)
# for str in unique:
#     print(variant_list.count(str))

