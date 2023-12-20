import gzip
import os
import re


class GtfFile():

    def __init__(self, gtf_path):
        self.gtf_path = gtf_path
        self.gencode_v = None
        self.genome_v = None
        self.ensembl_v = None
        self.date = None
        self.is_gz = gtf_path.endswith(".gz")
        self.header = self.get_header()
        self.get_header_versions()

    def open_file(self):
        if self.is_gz:
            return gzip.open(self.gtf_path, "rt")
        else:
            return open(self.gtf_path, "r")

    def get_header(self):
        """
        Get header from gtf file

        Returns:
        header(str): header of the gtf file
        """
        header = ""
        with self.open_file() as f:
            for line in f:
                if line.startswith("#"):
                    header += line
                # when non header line is found the function ends
                else:
                    return(header)
    
    def get_header_versions(self):
        if "GRCh37" in self.header:
            self.genome_v = 37
        elif "GRCh38" in self.header:
            self.genome_v = 38
        
        print(f"Gtf genome version: {self.genome_v}")
        
        # difining pattern to look for gencode version
        gencode_v_pattern = re.compile(r"version (\d+)")
        gencode_v_match = gencode_v_pattern.search(self.header)

        if gencode_v_match:
            self.gencode_v = gencode_v_match.group(1)
            print(f"Gtf gencode version = {self.gencode_v}")

        ensembl_v_pattern = re.compile(r"Ensembl (\d+)")
        ensembl_v_match = ensembl_v_pattern.search(self.header)

        if ensembl_v_match:
            self.ensembl_v = ensembl_v_match.group(1)
            print(f"Gtf Ensembl version = {self.ensembl_v}")
        
        # looking for string in format yyyy-mm-dd
        date_pattern = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")
        date_match = date_pattern.search(self.header)

        if date_match:
            self.date = date_match.group()
            print(f"Date of the gtf: {self.date}")


    def parse_gtf(self):
        """
        iterator, parses the gtf_path file and returns the line.strip()
        """
        with self.open_file() as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip()
                yield (line)
    
    def parse_gtf_line(self):
        """
        parses the gtf file and get different fields

        if gtf is sorted in location order it won't work as we need that for each gene, the sequence of elements found is:
        gene-transcrip-(exons/introns)-cds
        """
        # identify instance given an ensembl id
        genes_dict = dict()
        transcripts_dict = dict()
        exons_dict = dict()

        # sets where gtf objects will be stored
        genes_list = list()
        transcripts_list = list()
        exons_list = list()
        introns_list = list()
        cds_list = list()
        start_codon_list = list()
        stop_codon_list = list()
        utr_list = list()
        selenocysteine_list = list()


        genomic_type_set = set()

        for line in self.parse_gtf():
            line = line.split("\t")
            genomic_type = line[2]
            genomic_type_set.add(genomic_type)
            if genomic_type == "gene":
                # creating gene instance and added to genes dictionary
                gtf_gene_ins = GtfGene(line)
                genes_dict[gtf_gene_ins.gene_id] = gtf_gene_ins
                genes_list.append(gtf_gene_ins)

            elif genomic_type == "transcript":
                gtf_transcript_ins = GtfTranscript(line)
                if hasattr(gtf_transcript_ins, "tag"):
                    if gtf_transcript_ins.tag == "MANE_Select":
                        gtf_transcript_ins.is_mane = "MANE_Select"
                transcripts_dict[gtf_transcript_ins.transcript_id] = gtf_transcript_ins
                transcripts_list.append(gtf_transcript_ins)

            elif genomic_type == "exon":
                gtf_exon_ins = GtfExon(line)
                exon_id = gtf_exon_ins.exon_id.split("_")[0]
                exons_dict[exon_id] = gtf_exon_ins
                exons_list.append(gtf_exon_ins)

            elif genomic_type == "intron":
                gtf_intron_ins = GtfIntron(line)
                introns_list.append(gtf_intron_ins)

            elif genomic_type == "CDS":
                gtf_cds_ins = GtfCds(line)
                cds_list.append(gtf_cds_ins)

            elif genomic_type == "start_codon":
                gtf_start_codon_ins = GtfStartCodon(line)
                start_codon_list.append(gtf_start_codon_ins)

            elif genomic_type == "stop_codon":
                gtf_stop_codon_ins = GtfStopCodon(line)
                stop_codon_list.append(gtf_stop_codon_ins)

            elif genomic_type == "UTR":
                gtf_utr_ins = GtfUtr(line)
                utr_list.append(gtf_utr_ins)

            elif genomic_type == "Selenocysteine":
                gtf_selenocysteine_ins = GtfSelenocysteine(line)
                selenocysteine_list.append(gtf_selenocysteine_ins)
                # print(line)

        features_list = {
            "genes": genes_list,
            "transcripts": transcripts_list,
            "exons": exons_list,
            "introns": introns_list,
            "cds": cds_list,
            "start_codon": start_codon_list,
            "stop_codon": stop_codon_list,
            "utr": utr_list,
            "selenocysteine": selenocysteine_list
        }

        ids_gtf_object = {
            "genes": genes_dict,
            "transcripts": transcripts_dict,
            "exons": exons_dict,
        }
        return(features_list, ids_gtf_object)

    def connect_gtf_objs(self):
        features_list, ids_gtf_object = self.parse_gtf_line()

        genes_dict = ids_gtf_object["genes"]
        transcripts_dict = ids_gtf_object["transcripts"]
        exons_dict = ids_gtf_object["exons"]

        print("Associating each entity with its parent, e.g: gene->transcripts->exons/introns/CDS")

        for selenocysteine_ins in features_list["selenocysteine"]:
            transcript_id = selenocysteine_ins.transcript_id

            transcript_ins = transcripts_dict[transcript_id]
            transcript_ins.selenocysteine.append(selenocysteine_ins)
        
        for utr_ins in features_list["utr"]:
            transcript_id = utr_ins.transcript_id

            transcript_ins = transcripts_dict[transcript_id]
            transcript_ins.utr.append(utr_ins)

        
        for start_codon_ins in features_list["start_codon"]:
            exon_id = start_codon_ins.exon_id
            transcript_id = start_codon_ins.transcript_id

            exon_ins = exons_dict[exon_id]
            transcript_ins = transcripts_dict[transcript_id]

            exon_ins.start_codon.add(start_codon_ins)
            transcript_ins.start_codon.add(start_codon_ins)

        for stop_codon_ins in features_list["stop_codon"]:
            exon_id = stop_codon_ins.exon_id
            transcript_id = stop_codon_ins.transcript_id

            exon_ins = exons_dict[exon_id]
            transcript_ins = transcripts_dict[transcript_id]

            exon_ins.stop_codon.add(stop_codon_ins)
            transcript_ins.stop_codon.add(stop_codon_ins)
        
        for cds_ins in features_list["cds"]:
            exon_id = cds_ins.exon_id
            transcript_id = cds_ins.transcript_id

            exon_ins = exons_dict[exon_id]
            transcript_ins = transcripts_dict[transcript_id]

            exon_ins.cds.append(cds_ins)
            transcript_ins.cds.append(cds_ins)

        for intron_ins in features_list["introns"]:
            transcript_id = intron_ins.transcript_id
            gene_id = intron_ins.gene_id

            if transcript_id in transcripts_dict:
                transcript_ins = transcripts_dict[transcript_id]
            else:
                raise ValueError(
                    f"exon {exon_ins.exon_id} wich belongs to transcript id: {exon_ins.transcript_id}. The transcript ID "
                    "has not been found as a transcript object"
                )
            if gene_id in genes_dict:
                gene_ins = genes_dict[gene_id]
            else:
                raise ValueError(
                    f"gene_id {gene_id} not found in genes_dict"
                )

            transcript_ins.exons.append(intron_ins)
            gene_ins.add(intron_ins)

        for exon_ins in features_list["exons"]:
            gene_id = exon_ins.gene_id
            transcript_id = exon_ins.transcript_id

            if transcript_id in transcripts_dict:
                transcript_ins = transcripts_dict[transcript_id]
            else:
                raise ValueError(
                    f"exon {exon_ins.exon_id} wich belongs to transcript id: {exon_ins.transcript_id}. The transcript ID "
                    "has not been found as a transcript object"
                )
            if gene_id in genes_dict:
                gene_ins = genes_dict[gene_id]
            else:
                raise ValueError(
                    f"gene_id {gene_id} not found in genes_dict"
                )
            
            gene_ins.exons.add(exon_ins)
            transcript_ins.exons.append(exon_ins)
        
        # accessing each transcript object
        for transcript_ins in features_list["transcripts"]:
            gene_id = transcript_ins.gene_id
            gene_ins = genes_dict[gene_id]
            if transcript_ins.is_mane == "MANE_Select":
                gene_ins.mane_select_transcript = transcript_ins
            # if no transcript with most exons defined, set the current transcript
            if gene_ins.transcript_most_exons == None:
                gene_ins.transcript_most_exons = transcript_ins
            # if transcript with most exons already defined, 
            # compare the number of exons of the current transcript instance and the instance with most exons
            elif len(transcript_ins.exons) > len(gene_ins.transcript_most_exons.exons):
                gene_ins.transcript_most_exons = transcript_ins
            gene_ins.transcripts.append(transcript_ins)
    
        return(features_list)
        
    def create_spliceai_file(self):
        header = "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n"

        spliceai_dir = "/home/ocanal/ANN_DIR/spliceAI/hg19"
        filename = f"spliceAI_gencode{self.gencode_v}.txt"
        file_path = os.path.join(spliceai_dir, filename)
        print(file_path)

        features_list = self.connect_gtf_objs()
        genes = features_list["genes"]

        print(f"Creating new spliceAI file: {file_path}")
        with open(file_path, "w") as f:
            f.write(header)
            for gene in genes:
                if gene.mane_select_transcript:
                    transcript_ins = gene.mane_select_transcript
                else:
                    transcript_ins = gene.transcript_most_exons
                
                name = transcript_ins.gene_name
                chr = transcript_ins.chr

                if "chr" in chr:
                    chr = chr.replace("chr", "")
            
                strand = transcript_ins.strand
                start = transcript_ins.start
                end = transcript_ins.end
                exons = transcript_ins.exons
            
                exon_starts = list()
                # dict with key= exon start, value = exon end
                exon_ends_dict = dict()
            
                for exon in exons:
                    exon_starts.append(int(exon.start)-1)
                    exon_ends_dict[int(exon.start)-1] = exon.end
                if exon_starts:
                    # sorting exons by its start position
                    exon_starts.sort()
                    # craeting exons_end list sorted by exon start position
                    exon_ends = [exon_ends_dict[int(exon_start)] for exon_start in exon_starts]
                
                if len(exon_starts) != len(exon_ends):
                    raise(ValueError(
                        f"number of exons starts and ends are different: {len(exon_starts)}\t {len(exon_ends)}"
                    ))
            
                if exon_starts and exon_ends:
                    exon_start_string = ",".join(str(exon_start) for exon_start in exon_starts)
                    exon_end_string = ",".join(str(exon_end) for exon_end in exon_ends)
                    # SpliceAI example file always ends with ,
                    exon_start_string += ","
                    exon_end_string += ","
                else:
                    exon_start_string = ""
                    exon_end_string = ""
                
                
                line = [name, chr, strand, start, end, exon_start_string, exon_end_string]
                output_line = "\t".join(line)
                output_line += "\n"
                f.write(output_line)

        return file_path

    def sort_spliceAI_file(self, spliceAI_file_path):
        """
        Sorting the SpliceAI transcript file by chromosome and start position of the transcripts:
        
        Params:
            spliceAI_file_path(str): abs path where spliceAI file is.
        Returns:
            out_file_path(str): abs path where sorted spliceAI file is.
        """ 
        directory = os.path.dirname(spliceAI_file_path)
        filename = os.path.basename(spliceAI_file_path)
        out_filename = f"sorted_{filename}"
        out_file_path = os.path.join(directory, out_filename)

        print(f"sorting {spliceAI_file_path}")
        cmd = f"sort -k2,2 -k4,4n {spliceAI_file_path} > sorted_{out_file_path}"
        print(f"Sorting completed, new file in: {out_file_path}")
        subprocess.run(cmd)

            
class GtfLine():
    def __init__(self, line):
        self.chr = line[0]
        self.source = line[1]
        self.genomic_type = line[2]
        self.start = line[3]
        self.end =  line[4]
        self.score = line[5]
        self.strand = line[6]
        self.frame = line[7]
        self.attributes = line[8]
        self.parse_attributes()

    def parse_attributes(self):
        attributes = self.attributes.split(";")
        for attribute in attributes:
            if attribute == "":
                continue
            attribute = attribute.strip()
            key_value_attr = attribute.split(" ")
            key_attr = key_value_attr[0]
            value_attr = key_value_attr[1]
            # removing ""
            value_attr = value_attr.replace("\"", "")
            setattr(self, key_attr, value_attr)

class GtfGene(GtfLine):
    def __init__(self, line):
        super().__init__(line)
        self.transcripts = []
        self.exons = set()
        self.introns = set()
        self.mane_select_transcript = None
        self.transcript_most_exons = None

class GtfTranscript(GtfLine):
    def __init__(self, line):
        super().__init__(line)
        self.introns = []
        self.exons = []
        self.cds = []
        self.start_codon = set()
        self.stop_codon = set()
        self.selenocysteine = []
        self.utr = []
        self.is_mane = False

class GtfExon(GtfLine):
    def __init__(self, line):
        super().__init__(line)
        self.start_codon = set()
        self.stop_codon = set()
        self.cds = []
        self.utr = []

class GtfIntron(GtfLine):
    def __init__(self, line):
        super().__init__(line)

class GtfCds(GtfLine):          
    def __init__(self, line):
        super().__init__(line)
        self.selenocysteine = []

class GtfStartCodon(GtfLine):          
    def __init__(self, line):
        super().__init__(line)    

class GtfStopCodon(GtfLine):          
    def __init__(self, line):
        super().__init__(line)

class GtfUtr(GtfLine):          
    def __init__(self, line):
        super().__init__(line)

class GtfSelenocysteine(GtfLine):          
    def __init__(self, line):
        super().__init__(line)


gtf_ins = GtfFile("/home/ocanal/Desktop/parse_NGS_files/gencode.v44lift37.annotation.gtf.gz")
gtf_ins.create_spliceai_file()