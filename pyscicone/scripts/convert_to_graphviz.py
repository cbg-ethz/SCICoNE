import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("tree_file", help="The tree file output by SCICoNE.")
#parser.add_argument("--cell_node_ids_file", help="The attachment of each cell to a node", default="")
args = parser.parse_args()

tree_file = args.tree_file
#cell_node_ids_file = args.cell_node_ids_file

tree_filename = os.path.splitext(tree_file)[0]
output_file = tree_filename + '.graphviz'

root_label= "Neutral"
node_labels = True
color = "#E6E6FA"
event_fontsize=14
nodesize_fontsize=14
nodelabel_fontsize=14
node_sizes = None

with open(tree_file) as f:
	list_tree_file = list(f)

tree_str = ''.join(list_tree_file)

list_tree_file = tree_str.split('\n')

node_dict = dict()
for line in list_tree_file:
	if line.startswith("node"):
		colon_splits = line.split(':', 2)
		node_id = colon_splits[0].split(" ")[1]
		parent_id, event_str = colon_splits[2].split(',', 1)
		event_str = event_str[1:-1] # Remove '[' and ']'
		if parent_id != 'NULL':
			event_dict = dict([s.split(':') for s in event_str.split(',')])
		else:
			event_dict = dict()
		node_dict[node_id] = dict(parent_id=parent_id, region_event_dict=event_dict)



graphviz_header = [
		"digraph {",
		f'node [style=filled,color="{color}",fontsize={event_fontsize},margin=0,shape=oval]'
		f'edge [arrowhead=none, color="{color}"]',
	]

graphviz_labels = []
graphviz_links = []

for key in node_dict:
		node_id = key
		p_id = node_dict[key]['parent_id']
		if node_id == '0':
				str_merged_labels = root_label
		elif node_id == '-100':
				str_merged_labels = 'Whole-genome duplication'
		else:
				merged_labels = []

				event_dict = node_dict[key]['region_event_dict']
				first_region = list(event_dict.keys())[0]
				previous_event = event_dict[first_region]
				last_region = first_region
				for i, region in enumerate(event_dict):
						if i > 0:
							event = event_dict[region]
							if int(region) == int(last_region) + 1 and event == previous_event:
								last_region = region  # update the end
							else:
								if first_region == last_region:
									merged_labels.append(f"{int(previous_event):+}R{first_region}")
								else:
									merged_labels.append(f"{int(previous_event):+}R{first_region}:{last_region}")
								first_region = last_region = region
							previous_event = event
				if first_region == last_region:
					merged_labels.append(f"{int(previous_event):+}R{first_region}")
				else:
					merged_labels.append(f"{int(previous_event):+}R{first_region}:{last_region}")

				# Add line breaks
				region_str_merged_labels = " ".join(
						f"{x}<br/>" if i % 5 == 0 and i > 0 else str(x)
						for i, x in enumerate(merged_labels)
					)
				# Remove trailing line break
				if ''.join(list(region_str_merged_labels)[-len('<br/>'):]) == '<br/>':
					region_str_merged_labels = ''.join(list(region_str_merged_labels)[:-len('<br/>')])

				str_merged_labels = region_str_merged_labels

		# Add node size
		if node_sizes is not None:
			try:
				node_size = node_sizes[node_id]
			except KeyError:
				node_size = 0
			str_merged_labels = str_merged_labels + "<br/><br/>"
			str_merged_labels = (
			str_merged_labels
				+ f'<font point-size="{nodesize_fontsize}">'
				+ str(int(node_size))
				+ " cell"
			)
			if int(node_size) > 1 or int(node_size) == 0:
				str_merged_labels = str_merged_labels + "s"
			str_merged_labels = str_merged_labels + "</font>"

		graphviz_labels.append(
			f"{node_id}[label=<{str_merged_labels}>]"
		)  # use < > to allow HTML
		if p_id != 'NULL':
			graphviz_links.append(f"{p_id} -> {node_id}")

graphviz_str = '\n'.join(graphviz_header + graphviz_labels + graphviz_links + ["}"])

with open(output_file, 'w') as f:
	f.write(graphviz_str)
