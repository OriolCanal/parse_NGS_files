class Vcf_file():
    vep_info_line = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: "
    
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path

    def parse_vcf_header(self):
        with open(self.vcf_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    yield(line)

    def parse_vcf_lines(self):
        with open(self.vcf_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                fields = line.split("\t")

                yield(fields)
    
    def get_vep_header_annotations(self):

        for line in self.parse_vcf_header():
            if Vcf_file.vep_info_line in line:
                # removing description
                line = line.replace(Vcf_file.vep_info_line, "")
                # removing last chars ">
                line = line.replace("\">", "")
                # each vep field is seperated by |
                vep_desc_fields = line.split("|")

                return(vep_desc_fields)
    
    def get_isoforms(self):

        # obtaining the different fields that are annotated by vep
        vep_header_desc = self.get_vep_header_annotations()
        for line in self.parse_vcf_lines():
            info = line[7]
            info_fields = info.split(";")
            vep_field = next(field for field in info_fields if field.startswith("CSQ="))
            # removing CSQ= from the string
            vep_field = vep_field.replace("CSQ=", "")
            isoforms = vep_field.split(",")
            for isoform in isoforms:
                vep_desc_value = dict()
                isoform_vep_fields = isoform.split("|")
                for vep_desc, vep_value in zip(vep_header_desc, isoform_vep_fields):
                    vep_desc_value[vep_desc] = vep_value
                yield(vep_desc_value, line)

    def get_variant_isoforms(self):
        """
        get a dictionary with variants as keys in fomat chr:pos_ref>alt and the values are the vep 
        attrubutes associated with a variant
        """
        variants_isoforms = dict()
        for vep_desc_value, line in self.get_isoforms():

            chr = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]

            variant = f"{chr}:{pos}_{ref}>{alt}"
            if not variant in variants_isoforms:
                variants_isoforms[variant] = [vep_desc_value]
            else:
                variants_isoforms[variant].append(vep_desc_value)
        
        return(variants_isoforms)

    def compare_gff_vep_isoform(self):
        variants_isoforms_inf = self.get_variant_isoforms()
        no_common_fields = set()
        for variant, vep_desc_vep_result in variants_isoforms_inf.items(): 
            gencode_isofs = []
            vep_isofs = []
            for isoform in vep_desc_vep_result:
                if isoform["SOURCE"] == "GENCODE":
                    gencode_isofs.append(isoform)
                else:
                    vep_isofs.append(isoform)
            for gencode_iso in gencode_isofs:
                for vep_iso in  vep_isofs:
                    gencode_transid = gencode_iso["Feature"]
                    gencode_transid = gencode_transid.split(".")[0]
                    if gencode_transid != vep_iso["Feature"]:
                        continue
                    if gencode_iso["Allele"] != vep_iso["Allele"]:
                        continue
                    # print("___________----new isoform----_______")
                    for vep_desc in gencode_iso:
                        # if vep_desc != "REVEL_score":
                        #      continue
                        gencode_value = gencode_iso[vep_desc]
                        vep_value = vep_iso[vep_desc]

                        if not gencode_value == vep_value:
                            # print(variant)
                            # print(f"vep_desc: {vep_iso}\n\ngencode_desc: {gencode_iso}")
                            # print(vep_desc, gencode_value + " --- ", vep_value, "\n\n\n\n\n")
                            no_common_fields.add(vep_desc)
        return(no_common_fields)







    





rb36200_vcf = Vcf_file("/home/ocanal/Desktop/vcf_for_vep/RB36200.vep_gtf.vcf")
# for line in rb36200_vcf.parse_vcf_lines():
#     print(line)
vep_header = rb36200_vcf.get_vep_header_annotations()
# print(vep_header)
# for iso in rb36200_vcf.get_isoforms():

#     print(iso)
no_common_fields = rb36200_vcf.compare_gff_vep_isoform()
print(no_common_fields)


    

