import sys
import os
import re
import json
import hashlib

def hash_bytestr_iter(bytesiter, hasher, ashexstr=False):
    for block in bytesiter:
        hasher.update(block)
    return hasher.hexdigest() if ashexstr else hasher.digest()

def file_as_blockiter(afile, blocksize=65536):
    with afile:
        block = afile.read(blocksize)
        while len(block) > 0:
            yield block
            block = afile.read(blocksize)

def parse_config(config_fname):

    config = dict()
    for line in open(config_fname, 'r'):
        if line[0] == '#':
            continue
        if not '=' in line:
            continue
        key, val = line.strip().split('=')
        ### get rid of any comments
        val = re.sub(r'#.*$', '', val)
        ### replace bash variables (assume they have already been defined)
        for s in re.findall(r'\${.*\}', val):
            val = re.sub('\\' + s, config[s[2:-1]], val)
        ### handle bash arrays
        if val[0] == '(':
            val = [_.strip('"') for _ in val[1:-1].split(' ')]
        config[key] = val
    return config

def main():
    
    if len(sys.argv) < 2:
        sys.stderr.write('Usage: %s <config>\n' % sys.argv[0])
        sys.exit(1)
    config_fname = sys.argv[1]
    config = parse_config(config_fname)

    ### make outdir
    outdir = os.path.join(config['BASEDIR'], config['SIM_ID'], 'json')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    payload = json.load(open('templates/template.payload.json', 'r'))
    job_star = json.load(open('templates/STAR.template.job.json', 'r'))
    job_hisat2 = json.load(open('templates/HISAT2.template.job.json', 'r'))

    payload['studyId'] = config['SIM_ID']
    for i, genome in enumerate(config['GENOMES']):
        sid = genome + '_' + config['ICGC_SAMPLES'][i]
        for tn in ['tumor', 'normal']:
            for n in range(1, int(config['NUM_SAMPLES']) + 1):
                payload['experiment']['submitter_sequencing_experiment_id'] = sid
                payload['read_groups'][0]['file_r1'] = os.path.join(config['BASEDIR'], config['SIM_ID'], 'reads', f'{sid}_{tn}_sample_0{n}_1.fastq.gz')
                payload['read_groups'][0]['file_r2'] = os.path.join(config['BASEDIR'], config['SIM_ID'], 'reads', f'{sid}_{tn}_sample_0{n}_2.fastq.gz')
                payload['read_groups'][0]['library_name'] = f'{sid}_{tn}_sample_0{n}'
                payload['samples'][0]['sampleId'] = f'{sid}_{tn}_sample_0{n}'
                payload['samples'][0]['specimenId'] = f'{sid}_{tn}_sample_0{n}'
                payload['samples'][0]['submitterSampleId'] = f'{sid}_{tn}_sample_0{n}'
                payload['samples'][0]['specimen']['specimenId'] = f'{sid}_{tn}'
                payload['samples'][0]['specimen']['donorId'] = f'{sid}'
                payload['samples'][0]['specimen']['submitterSpecimenId'] = f'{sid}_{tn}'
                payload['samples'][0]['specimen']['tumourNormalDesignation'] = tn
                payload['samples'][0]['donor']['donorId'] = sid
                payload['samples'][0]['donor']['submitterDonorId'] = sid
                payload['files'][0]['fileName'] = os.path.join(config['BASEDIR'], config['SIM_ID'], 'reads', f'{sid}_{tn}_sample_0{n}_1.fastq.gz')
                payload['files'][0]['fileSize'] = os.path.getsize(payload['files'][0]['fileName'])
                payload['files'][0]['fileMd5sum'] = hash_bytestr_iter(file_as_blockiter(open(payload['files'][0]['fileName'], 'rb')), hashlib.md5(), ashexstr=True) 
                payload['files'][1]['fileName'] = os.path.join(config['BASEDIR'], config['SIM_ID'], 'reads', f'{sid}_{tn}_sample_0{n}_2.fastq.gz')
                payload['files'][1]['fileSize'] = os.path.getsize(payload['files'][1]['fileName'])
                payload['files'][1]['fileMd5sum'] = hash_bytestr_iter(file_as_blockiter(open(payload['files'][1]['fileName'], 'rb')), hashlib.md5(), ashexstr=True) 

                with open(os.path.join(outdir, f'{sid}_{tn}_sample_0{n}.payload.json'), 'w', encoding='utf-8') as out:
                    json.dump(payload, out, ensure_ascii=False, indent=4)

                job_star['index'] = config['INDEX_STAR']
                job_star['gtf'] = config['ANNO_GTF']
                job_star['input_files'][0] = payload['read_groups'][0]['file_r1']
                job_star['input_files'][1] = payload['read_groups'][0]['file_r2']
                job_star['metadata'] = os.path.join(outdir, f'{sid}_{tn}_sample_0{n}.payload.json')
                job_star['publish_dir'] = config['ALIGN_OUTDIR'] 
                job_star['container'] = config['STAR_IMAGE']

                with open(os.path.join(outdir, f'{sid}_{tn}_sample_0{n}.STAR.job.json'), 'w', encoding='utf-8') as out:
                    json.dump(job_star, out, ensure_ascii=False, indent=4)

                job_hisat2['index'] = config['INDEX_HISAT2']
                job_hisat2['gtf'] = config['ANNO_GTF']
                job_hisat2['input_files'][0] = payload['read_groups'][0]['file_r1']
                job_hisat2['input_files'][1] = payload['read_groups'][0]['file_r2']
                job_hisat2['metadata'] = os.path.join(outdir, f'{sid}_{tn}_sample_0{n}.payload.json')
                job_hisat2['publish_dir'] = config['ALIGN_OUTDIR'] 
                job_hisat2['container'] = config['HISAT2_IMAGE']

                with open(os.path.join(outdir, f'{sid}_{tn}_sample_0{n}.HISAT2.job.json'), 'w', encoding='utf-8') as out:
                    json.dump(job_hisat2, out, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()

