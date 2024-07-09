import json

def output_acc_org_asm(filename):
    with open(filename, 'r') as handle:
        for i, line in enumerate(handle, 1):
            record = json.loads(line)
            acc = record['accession']
            org = record['organism']['organismName'].replace(' ', '_')
            asm = record['assemblyInfo']['assemblyName'].replace(' ', '_')
            print(i, acc, org, asm)