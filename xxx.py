from pathlib import Path


def remove_with_name(s: str, name: str) -> str:
    columns = s.split('\t')
    last_column = columns[-1]
    last_columns = last_column.split(';')
    last_columns = [x for x in last_columns if not x.startswith(name + '=')]
    last_column = ';'.join(last_columns)
    columns[-1] = last_column
    return '\t'.join(columns)


def generate(vcf_path: str):
    names = ['DP', 'SB', 'DP4', 'AF']
    new_path = vcf_path.replace('.vcf', f'_no_info.vcf')
    print('output to ->', new_path)
    with open(new_path, 'w') as f:
        for line in Path(vcf_path).read_text().splitlines():
            if line.startswith('#'):
                f.write(line + '\n')
                continue
            line = line.strip()
            for name in names:
                line = remove_with_name(line, name)
            f.write(line + '\n')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', type=str, required=True, help='vcf path')
    args = parser.parse_args()
    vcf_path = args.path
    generate(vcf_path)
