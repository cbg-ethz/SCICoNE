import string
import seaborn as sns
from matplotlib.colors import ListedColormap

LABEL_CPAL = sns.color_palette('Set2', 20)
LABEL_CPAL_HEX = LABEL_CPAL.as_hex()
LABEL_CMAP = ListedColormap(LABEL_CPAL_HEX)
LABEL_COLORS_DICT = dict(zip(list(string.ascii_uppercase)[:20], LABEL_CPAL_HEX))

BLUE_WHITE_RED = ["#2040C8", "white", "#EE241D"]
