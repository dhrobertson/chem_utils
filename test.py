import pprint
from chembl_webresource_client.new_client import new_client
target = new_client.target
gene_name = 'SHIP1'
res = target.search(gene_name)
print(len(res))
print(pprint.pformat(res))
