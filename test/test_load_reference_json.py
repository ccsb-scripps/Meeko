from meeko import Polymer
fn = "/Users/amyhe/Desktop/0_forks/Meeko/test/polymer_data/AHHY_reference_fewer_templates.json"
with open(fn) as f:
	json_string = f.read()
polymer = Polymer.from_json(json_string)
polymer.get_valid_monomers()
