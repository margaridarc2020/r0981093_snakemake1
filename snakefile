DB = "/staging/leuven/stg_00079/teaching/hg38_21/chr21.fa"
SNPEFF_JAR="/lustre1/project/stg_00079/teaching/I0U19a_conda_2025/share/snpeff-5.2-1/snpEff.jar"

import glob
all_samples = glob_wildcards('000.fastq/{sample}.fastq').sample

rule all:
    default_target: True
    input:
        expand('010.fastqc/{sample}_fastqc.zip', sample=all_samples),
        ['020.bwa/{}_bwa.bam'.format(sample) for sample in all_samples],
        expand('020.bwa/{sample}_bwa.bam.bai',sample=all_samples),
        expand('030.samtools/{sample}_stats.tsv', sample=all_samples),        
        '030.samtools/raw_snps.vcf',
        '040.cleaned/clean_snps.vcf',
        '050.snpeff/annotated_snps.vcf',
        'vcf.sqlite',
        'high_impact_differential_snps.tsv'

rule fastqc:
	input:
		'000.fastq/{sample}.fastq'
	output:
		'010.fastqc/{sample}_fastqc.zip',
	shell:
		'''
        mkdir -p 010.fastqc &&
        fastqc {input} -o 010.fastqc/ 
        '''

rule bwa:
    input:
        '000.fastq/{sample}.fastq'
    output:
        '020.bwa/{sample}_bwa.bam'
    shell:
        '''
        mkdir -p 020.bwa &&
        bwa mem {DB} {input} | samtools sort -o {output} -        
        '''

rule bwa_bai:
    input:
        '020.bwa/{sample}_bwa.bam'
    output:
        '020.bwa/{sample}_bwa.bam.bai'
    shell:
        'samtools index {input} {output}'

rule samtools_stats:
    input:
        '020.bwa/{sample}_bwa.bam'
    output:
        '030.samtools/{sample}_stats.tsv'
    shell:
        '''
        mkdir -p 030.samtools &&
        samtools flagstat {input} > {output}
        '''

rule samtools:
    input:
        '020.bwa/TLE66_N_bwa.bam', 
        '020.bwa/TLE66_T_bwa.bam'
    output:
        '030.samtools/raw_snps.vcf'
    shell:
        '''
        mkdir -p 030.samtools &&
        bcftools mpileup -Ou -f {DB} {input} | bcftools call -mv -Ov -o {output}
        '''

rule cleaned:
    input:
        '030.samtools/raw_snps.vcf'
    output:
        '040.cleaned/clean_snps.vcf'
    shell:
        '''
        mkdir -p 040.cleaned &&
        cat {input} | vt decompose - \
        | vt normalize -n -r {DB} - \
        | vt uniq - | vt view -f "QUAL>20" -h - > {output}
        '''

rule snpeff:
    input:
        '040.cleaned/clean_snps.vcf'
    output:
        '050.snpeff/annotated_snps.vcf'
    shell:
        '''
        mkdir -p 050.snpeff &&
        java -Xmx3400m -jar {SNPEFF_JAR} eff hg38 \
        -dataDir /staging/leuven/stg_00079/teaching/snpeff_db \
        {input} \
        > {output}
        '''

rule vcf_db:
    input:
        vcf="050.snpeff/annotated_snps.vcf"
    output:
        db="vcf.sqlite"
    run:
        import pandas as pd
        import vcfpy
        import os 
        import sqlite3
        import sys
        
        vcf = input.vcf
        db = output.db
        
        dbfile = os.path.join(db)
        db = sqlite3.connect(dbfile)      
        
        snp_records = []
        call_records = []
        effect_records = []
        
        reader = vcfpy.Reader.from_path(vcf)
        
        effect_rec_names = """snp allele effect impact gene gene_id feature_type feature_id biotype rank hgvs.c hgvs.p cdna_pos cds_pos 
                              prot_pos distance_to_feature messages""".split()
        
        for i, record in enumerate(reader):
            assert len(record.ALT) == 1
            alt = record.ALT[0]
        
            snp_name = f"{record.CHROM}:{record.POS}:{record.REF}:{alt.value}"
        
            snp_records.append(
                dict(snp = snp_name,
                     chrom = record.CHROM,
                     pos = record.POS,
                     qual = record.QUAL,
                     ref = record.REF,
                     type = alt.type,
                     alt = alt.value))
        
            for call_record in record.calls:
                sample = os.path.basename(call_record.sample)
                sample = sample.replace('.bam', '')
                gt = call_record.data['GT']
                if gt == '1/0':
                    gt = '0/1'
        
                call_records.append(
                    dict(snp=snp_name,
                         sample = sample, 
                         patient = sample.split('_')[0],
                         status = sample.split('_')[1],
                         genotype = gt))
        
            for ann in record.INFO['ANN']:
                ann = ann.split('|')
               
                eff_record = dict(zip(effect_rec_names, [snp_name] + ann))
                try: 
                    eff_record['distance_to_feature'] = int(eff_record.get(intfield)) 
                except: 
                    eff_record['distance_to_feature'] = -1 
                
                effect_records.append(eff_record)
        
        snp_records = pd.DataFrame.from_records(snp_records)
        call_records = pd.DataFrame(call_records)
        effect_records = pd.DataFrame(effect_records)
        
        print('snp records:', snp_records.to_sql('snp', db, if_exists = 'replace', index=False))  
        print('snp calls:', call_records.to_sql('snp_call', db, if_exists='replace', index=False))
        print('snp effects:', effect_records.to_sql('snp_effect', db, if_exists='replace', index=False))
        
        db.commit()
        db.close()


rule sqlite_to_tsv:
    input:
        vcf="vcf.sqlite"
    output:
        tsv="high_impact_differential_snps.tsv"
    run:
        import pandas as pd
        import sqlite3
        import sys

        db_path = input.vcf
        tsv_path = output.tsv
        
        db = sqlite3.connect(db_path)
        
        read_call_effect = pd.read_sql("""select distinct snp.snp, snp.chrom, snp.pos, snp.alt, call_N.genotype as ng, call_T.genotype as tg, snp_effect.impact
            from snp, snp_call, snp_call as call_N, snp_call as call_T, snp_effect
            where snp.snp = call_N.snp
                and snp.snp = call_T.snp
                and call_N.status = 'N'
                and call_T.status = 'T'
                and call_N.genotype != call_T.genotype
                and NOT call_N.genotype LIKE '%.%'
                and NOT call_T.genotype LIKE '%.%'
                and snp_effect.impact == 'HIGH'
                and snp_call.snp = snp_effect.snp
            """, db)
        
        read_call_effect.to_csv(tsv_path, sep='\t', index=False)
        db.close()


