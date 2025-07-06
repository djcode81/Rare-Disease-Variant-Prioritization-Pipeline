#!/usr/bin/env python3

import vcf

def debug_vep_parsing(vcf_path):
    print(f"🔍 Debugging VEP parsing for: {vcf_path}")
    print("=" * 60)
    
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    
    for i, record in enumerate(vcf_reader):
        print(f"\n📍 Record {i+1}:")
        print(f"   Position: {record.CHROM}:{record.POS}")
        print(f"   REF>ALT: {record.REF}>{record.ALT}")
        
        # Check if CSQ exists
        if 'CSQ' in record.INFO:
            csq_data = record.INFO['CSQ']
            print(f"   CSQ type: {type(csq_data)}")
            print(f"   CSQ data: {csq_data}")
            
            if isinstance(csq_data, list):
                print(f"   CSQ list length: {len(csq_data)}")
                if len(csq_data) > 0:
                    first_csq = csq_data[0]
                    print(f"   First CSQ: {first_csq}")
                    
                    fields = first_csq.split('|')
                    print(f"   Fields count: {len(fields)}")
                    print(f"   Field 0 (ALT): {fields[0] if len(fields) > 0 else 'MISSING'}")
                    print(f"   Field 1 (Consequence): {fields[1] if len(fields) > 1 else 'MISSING'}")
                    print(f"   Field 2 (Impact): {fields[2] if len(fields) > 2 else 'MISSING'}")
                    print(f"   Field 3 (SYMBOL): {fields[3] if len(fields) > 3 else 'MISSING'}")
                    
                    # Look for scores in later fields
                    for j, field in enumerate(fields):
                        if field and ('.' not in field or field.replace('.', '').isdigit()):
                            if j > 15:  # Scores usually later in VEP format
                                print(f"   Field {j}: {field}")
            else:
                print(f"   CSQ is not a list: {csq_data}")
        else:
            print("   ❌ No CSQ field found")
    
    print("\n" + "=" * 60)

if __name__ == "__main__":
    debug_vep_parsing('test_data/real_vep_annotated.vcf')
