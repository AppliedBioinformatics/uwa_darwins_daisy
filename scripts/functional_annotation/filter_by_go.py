import pandas as pd
import plotly.graph_objects as go

# Function to extract the uniprot accession from "subject_id".
def extract_accession(s):
    # If format is like "sp|F8S1H3|C7BL1_HELAN", extract middle part
    if '|' in s:
        parts = s.split('|')
        if len(parts) > 1:
            return parts[1]
    # Otherwise return the string as-is
    return s

# Import Non-core Go Merged gene list.
df = pd.read_csv("../../data/functional_annotation/noncore_go_merged_diamond_results_uniprot.csv")

# GO TERMS TO FILTER.
go_terms = ["GO:0016114"]
filtered_df = df[df["go_terms"].str.contains("GO:0016114", na=False)]

# Generate summary table.
summary = filtered_df.groupby("subject_id").agg(
    n_copies_in_gff=("subject_id", "count"),
    avg_bit_score=("bit_score", "mean"),
    reference_species=('sscinames', lambda x: ', '.join(sorted(set(x)))),
    description=('stitle', "first"),
).reset_index().sort_values("n_copies_in_gff", ascending=False)

# Add uniprot links.
summary['uniprot_accession'] = summary['subject_id'].apply(extract_accession)
summary['uniprot_link'] = summary['uniprot_accession'].apply(
    lambda x: f'https://www.uniprot.org/uniprot/{x}'
)

# Convert DataFrame to HTML (escape=False to keep the links active)
html_table = summary.to_html(index=False, escape=True)

# DataTables CSS + JS CDN links
datatables_css = '<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css"/>'
datatables_js = '<script type="text/javascript" src="https://code.jquery.com/jquery-3.5.1.js"></script>' \
                '<script type="text/javascript" src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>'

# JavaScript to activate DataTables on our table
datatables_init = """
<script type="text/javascript">
$(document).ready(function() {
    $('table').DataTable({
        "paging": true,
        "searching": true,
        "info": true
    });
});
</script>
"""

# Full HTML page combining all parts
html_page = f"""
<html>
<head>
    <meta charset="utf-8" />
    {datatables_css}
</head>
<body>
    {html_table}
    {datatables_js}
    {datatables_init}
</body>
</html>
"""

# Save to file
with open("../../plots/functional_annotation/non_core_GO0016114_summary_table.html", "w", encoding="utf-8") as f:
    f.write(html_page)

