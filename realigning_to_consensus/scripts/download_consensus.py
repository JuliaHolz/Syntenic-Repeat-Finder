import requests
import argparse 


parser = argparse.ArgumentParser(
                    prog='download_consensus',
                    description='downloads consensus from Dfam',
                    epilog='using Dfam REST API')
parser.add_argument('-f', '--family', action="store", help = "family name")


args = parser.parse_args()

family_name = args.family

url = "https://dfam.org/api/families"
params = {
    # The summary format is metadata-only and does not include
    "name":family_name,
    "format": "full",

}
response = requests.get(url, params)
results = response.json()["results"]
print(results[0].keys())
print(results[0]['consensus_sequence'])
family_consensus = results[0]['consensus_sequence']

consensus_file = "./family_consensi/"+family_name + ".fasta"
with open(consensus_file, "w") as cf:
    cf.write(">" + family_name + "_consensus\n")
    cf.write(family_consensus)