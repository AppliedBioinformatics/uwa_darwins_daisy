pop_colors = {
    'HAMP28_35': '#bbbc40',
    'HALM35_20': '#674ea7',
    'HALM36_12': '#cc0000',
    'HALM51_21': '#c27ba0',
    'HALM44_17': '#c27ba0',
    'HALM53_3':  '#c27ba0',
    'HAMP29_50': '#bbbc40',
    'HAMP10_45': '#6aa84f',
    'HAMP4_56':  '#6aa84f',
    'HALM62_6':  '#6aa84f',
    'HALM41_26': '#6aa84f',
    'HALM65_11': '#999999',
    'HALM26_18': '#6aa84f',
    'LAM26_15':  '#6aa84f',
    'HALM64_10': '#6aa84f',
    'HALM37_1':  '#660000',
    'HALM8_13':  '#6aa84f',
    'HALM39_23': '#6aa84f',
    'HAMP7_6':   '#6aa84f',
    'HAMP5_15':  '#6aa84f',
    'HALM19_24': '#f69035',
    'HALM11_13': '#f69035',
    'HAMP40_29': '#f69035',
    'HALM12_19': '#f69035',
    'HALM23_8':  '#f69035',
    'HALM20_21': '#f69035',
    'AHA6_30':   '#0b5394',
    'SLAM20_2':  '#9fc5e8',
    'SLAM27_24': '#9fc5e8',
    'HAMP38_7':  '#f69035',
    'HALM32_22': '#cc0000',
    'HALM25_19': '#6aa84f',
    'HAMP30_39': '#bbbc40',
    'HAMP36_1':  '#f69035',
}

order = [
    "HALM35_20", "HALM36_12", "HALM32_22", "HALM44_17", "HALM51_21", "HALM53_3", "AHA6_30",
    "HALM37_1", "HALM26_18", "HALM62_6", "HALM64_10",
    "HAMP10_45", "HAMP4_56", "HALM8_13", "HAMP5_15", "HAMP7_6", "HALM39_23", "HALM41_26", "LAM26_15", "HALM25_19", "HALM65_11",
    "SLAM20_2", "SLAM27_24", "HALM11_13", "HALM12_19", "HALM19_24", "HALM20_21", "HALM23_8",
    "HAMP40_29", "HAMP36_1", "HAMP38_7", "HAMP28_35", "HAMP29_50", "HAMP30_39"
]

ordered_pop_colors = {k: pop_colors[k] for k in order}

island_colors = pop_colors = {
    "Bartolomé": "#674ea7",
    "Santiago": "#cc0000",
    "Isabela": "#c27ba0",
    "Pinzon": "#0b5394",
    "Isla Eden": "#660000",
    "Santa Cruz": "#6aa84f",
    "Santa Fé": "#999999",
    "Pinta": "#9fc5e8",
    "San Cristóbal": "#f69035",
    "Floreana": "#bbbc40"
}