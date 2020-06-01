import string
import seaborn as sns

LABEL_CMAP = sns.color_palette('Set2', 20).as_hex()
LABEL_COLORS_DICT = dict(zip(list(string.ascii_uppercase)[:20], LABEL_CMAP))

BLUE_WHITE_RED = ["#2040C8", "white", "#EE241D"]
