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

    def parse_vcf_lines(self, header = False):
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
    
if __name__ == "__main__":
    rb36200_vcf = Vcf_file("/home/ocanal/Desktop/vcf_for_vep/RB36200.vep_gtf.vcf")

    vep_header = rb36200_vcf.get_vep_header_annotations()

    print(vep_header)




    

