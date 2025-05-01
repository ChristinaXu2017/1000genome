process CREATE_PCA {
    input:
    path vcf_file
    path panel_file
    output:
    path "pca_pop_plot.png"
    path "pca_pop_plot.html"
    path "pca_result.eigenvec"
    script:
    """
    plink --vcf ${vcf_file} --make-bed --out plink_data
    plink --bfile plink_data --pca 2 --out pca_result
    python3 - <<EOF
    import pandas as pd
    import plotly.express as px
    pca_df = pd.read_csv("pca_result.eigenvec", sep="\\s+", header=None, names=["FID", "IID", "PC1", "PC2"])
    panel_df = pd.read_csv("${panel_file}", sep='\\t')
    plot_df = pca_df.merge(panel_df[['sample', 'pop', 'super_pop']], left_on="IID", right_on="sample", how='left')
    fig = px.scatter(plot_df, x='PC1', y='PC2', title="PCA with Super Population Info", color='super_pop')
    fig.update_layout(height=800, width=1200,
                      xaxis=dict(showgrid=False, zeroline=False),
                      yaxis=dict(showgrid=False, zeroline=False))
    fig.write_image("pca_pop_plot.png")
    fig.write_html("pca_pop_plot.html")
    EOF
    """
}
