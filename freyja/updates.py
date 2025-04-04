import urllib.request
import os
import sys
import subprocess
import requests


def download_tree(locDir):
    # url = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/"\
    #       "UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
    # url = "https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/"\
    #       "UShER_SARS-CoV-2/2023/04/10/public-2023-04-10.all.masked.pb.gz"
    # url = "https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/"\
    #       "UShER_SARS-CoV-2/2023/08/17/public-2023-08-17.all.masked.nextclade.pangolin.pb.gz"
    #url = "https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/"\
    #        "UShER_SARS-CoV-2/2021/11/15/"\
    #        "public-2021-11-15.all.masked.nextclade.pangolin.pb.gz"
    #url = "https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/"\
    #        "UShER_SARS-CoV-2/2022/05/15/public-2022-05-15.all.masked.pb.gz"
    url = "https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/"\
            "UShER_SARS-CoV-2/2024/11/04/public-2024-11-04.all.masked.pb.gz"
    treePath = os.path.join(locDir, "public-latest.all.masked.pb.gz")
    urllib.request.urlretrieve(url, treePath)
    return treePath


def download_barcodes_wgisaid(locDir):
    url = "https://raw.githubusercontent.com/andersen-lab/"\
          "Freyja/main/freyja/data/usher_barcodes_with_gisaid.csv"
    bPath = os.path.join(locDir, "usher_barcodes_with_gisaid.csv")
    urllib.request.urlretrieve(url, bPath)
    return bPath


def download_barcodes(locDir):
    url = "https://raw.githubusercontent.com/andersen-lab/"\
          "Freyja/main/freyja/data/usher_barcodes.csv"
    bPath = os.path.join(locDir, "usher_barcodes.csv")
    urllib.request.urlretrieve(url, bPath)
    url2 = "https://raw.githubusercontent.com/andersen-lab/"\
           "Freyja/main/freyja/data/last_barcode_update.txt"
    bPath = os.path.join(locDir, "last_barcode_update.txt")
    urllib.request.urlretrieve(url2, bPath)
    return bPath


def convert_tree(loc_dir):
    print(f"Writing updated files to: {loc_dir}")
    tree_path = os.path.join(loc_dir, "public-latest.all.masked.pb.gz")
    var_cmd = f"matUtils extract -i {tree_path} -C lineagePaths.txt"
    sys.stdout.flush()  # force python to flush
    return_code = subprocess.run(var_cmd, shell=True, executable="/bin/bash",
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.PIPE)
    return return_code


def get_curated_lineage_data(locDir):
    # url2 = "https://raw.githubusercontent.com/outbreak-info/outbreak.info/"\
    #    "master/web/src/assets/genomics/curated_lineages.json"
    url2 = "https://raw.githubusercontent.com/outbreak-info/outbreak.info/"\
        "e6777f8c190f8f328b4450d6f9f8d836172cc82b/web/src/assets/genomics/"\
        "curated_lineages.json"
    urllib.request.urlretrieve(url2,
                               os.path.join(locDir,
                                            "curated_lineages.json"))


def get_cl_lineages(locDir):
    # for now, use lineages metadata created using patch
    r = requests.get('https://raw.githubusercontent.com/outbreak-info/' +
                     'outbreak.info/master/curated_reports_prep/lineages.yml')
    # r = requests.get('https://raw.githubusercontent.com/cov-lineages' +
    #                  '/lineages-website/master/data/lineages.yml')
    if r.status_code == 200:
        with open(os.path.join(locDir, 'lineages.yml'), 'w+') as f:
            f.write(r.text)


if __name__ == '__main__':
    download_tree()
    get_curated_lineage_data()
    convert_tree()
