import streamlit as st
import pandas as pd
import obonet
import requests
import networkx as nx
import re

# st.sidebar.subheader("OMIM diseases and associated HPO terms!")


def gene():
    st.title("Gene to Phenotype associations")
    # st.write(
    #    "Ditto+Hazel together assist clinical analysts and MDx lab directors in their rare disease genomic variant interpretation workflow by providing a ranked list of genes based on the patient's clinically recorded phenotype and by determining when there are likely to be more variants not yet identified that are contributing to the patients phenotype."
    # )

    @st.cache(allow_output_mutation=True)
    def load_data():
        # get HPO terms for OMIM disorders and write a output file mapping OMIM IDs to HPO terms
        graph = obonet.read_obo("http://purl.obolibrary.org/obo/hp.obo")
        inv_graph = graph.reverse()
        # Mapping from term ID to name
        id_to_name = {id_: data.get("name") for id_, data in graph.nodes(data=True)}

        res = {key + " : " + id_to_name[key]: key for key in id_to_name.keys()}

        gene2hpo = {}
        gene2hpo_assoc = {}
        r = requests.get(
            "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt", stream=True
        )
        # r.text
        for line in r.iter_lines():
            line = line.decode("utf-8")
            if not line.startswith("#"):
                part = line.strip("\n").split("\t")
                hp_term = part[2]
                # hp_term = alt_id.get(hp_term, hp_term)
                gene_name = part[1].upper()
                if gene_name not in gene2hpo:
                    gene2hpo[gene_name] = list([])
                    gene2hpo_assoc[gene_name] = list([])
                gene2hpo[gene_name].append(hp_term)
                gene2hpo_assoc[gene_name] = list(
                    set(gene2hpo_assoc[gene_name] + list(inv_graph.predecessors(hp_term)))
                    - set(gene2hpo[gene_name])
                )

                gene2hpo[gene_name] = list(
                    dict.fromkeys(gene2hpo[gene_name])
                )  # check for duplicate values and drop them
                gene2hpo[gene_name] = list(filter(None, gene2hpo[gene_name]))

                gene2hpo_assoc[gene_name] = list(
                    dict.fromkeys(gene2hpo_assoc[gene_name])
                )  # check for duplicate values and drop them
                gene2hpo_assoc[gene_name] = list(filter(None, gene2hpo_assoc[gene_name]))

        # for gene_name in gene2hpo.keys():
        #    gene2hpo_assoc[gene_name] = gene2hpo_assoc[gene_name] - gene2hpo[gene_name]

        gene2hpo_assoc_dict = {
            "Genes": list(gene2hpo_assoc.keys()),
            "Associated HPOs": list(gene2hpo_assoc.values()),
        }
        gene2hpo_assoc_dict = pd.DataFrame(gene2hpo_assoc_dict)

        gene2hpo_dict = {"Genes": list(gene2hpo.keys()), "HPOs": list(gene2hpo.values())}
        gene2hpo_dict = pd.DataFrame(gene2hpo_dict)

        gene2hpo_df = gene2hpo_dict.merge(gene2hpo_assoc_dict, on="Genes")
        gene2hpo_df["all"] = gene2hpo_df["HPOs"] + gene2hpo_df["Associated HPOs"]

        del (
            gene2hpo_assoc,
            gene2hpo,
            gene2hpo_dict,
            gene2hpo_assoc_dict,
            graph,
            inv_graph,
            part,
            hp_term,
            gene_name,
        )

        # gene2hpo_dict = gene_df(gene2hpo,graph)
        return gene2hpo_df, id_to_name, res  # ,gene2hpo_assoc

    @st.cache(allow_output_mutation=True)
    def gene_ranks(hpo_terms, gene_scores):
        term_len = len(hpo_terms)

        gene_scores["score"] = [
            (len(list(set(hpo_terms) & set(i)))) / term_len for i in gene_scores["all"]
        ]
        gene_scores = (
            gene_scores[["Genes", "score"]]
            .sort_values(by="score", ascending=False)
            .reset_index(drop=True)
        )

        # gene_scores.drop(["HPOs","Associated HPOs"], axis=1, inplace=True)
        gene_scores.index = gene_scores.index + 1
        return gene_scores

    data_load_state = st.text("Loading data...")
    # gene2hpo, gene2hpo_dict, graph, id_to_name, res = load_data()
    gene2hpo_df, id_to_name, res = load_data()
    data_load_state.text("Loaded HPO data!")
    data_load_state.text("")

    col1, col2, col3 = st.columns(3)

    # if col1.checkbox("Show loaded HPO metadata data"):
    #    col1.write(data.graph)
    terms = []
    options = col1.multiselect("Single HPO terms input use:", list(res.keys()))
    terms = [res[term] for term in options]
    hpo_text = col1.text_input("Input list of HPO terms seperated by comma/space/semi-colon", "")
    terms = terms + re.split(r"[,;\s]", hpo_text)
    terms = list(filter(None, terms))

    gene = col3.selectbox("Select a Gene:", gene2hpo_df.Genes.values)
    col2.write("Hazel scores for Genes:")

    # if col1.button('Submit') and len(terms) != 0:
    if len(terms) != 0:
        term_len = len(terms)
        list_terms = {key + " : " + id_to_name[key] for key in terms}
        col1.write(f"Total HPO terms = {term_len}")
        col1.write(list(list_terms))

        gene_scores = gene_ranks(terms, gene2hpo_df)
        cutoff = col2.slider('Score cutoff?', 0.0, 1.0, 0.1)
        gene_scores = gene_scores[gene_scores.score>cutoff]
        col2.download_button( "Download Gene ranks", gene_scores.to_csv(), "Hazel_predictions.csv", "text/csv", key="download-csv")
        gene_scores = gene_scores.style.format(subset="score", precision=2).bar(
            subset="score", align="mid"
        )
        col2.write(gene_scores.to_html(), unsafe_allow_html=True)

        gene_direct_hpo = list(
            set(terms) & set(gene2hpo_df["HPOs"].loc[gene2hpo_df["Genes"] == gene].to_list()[0])
        )
        gene_indirect_hpo = list(
            set(terms)
            & set(gene2hpo_df["Associated HPOs"].loc[gene2hpo_df["Genes"] == gene].to_list()[0])
        )
        score = (len(list(set(terms) & set(gene_direct_hpo + gene_indirect_hpo)))) / term_len
        col3.subheader(f"Hazel score = {score}")

        direct = {key + " : " + id_to_name[key] for key in gene_direct_hpo}
        indirect = {key + " : " + id_to_name[key] for key in gene_indirect_hpo}
        # all_associations = {
        #    key + " : " + id_to_name[key]
        #    for key in list(set(terms) & set(gene2hpo[gene] + gene2hpo_assoc))
        # }

        col3.write(f"\n\nDirect associations: {len(direct)}")
        col3.write(list(direct))
        col3.write(f"\n\nIndirect associations: {len(indirect)}")
        col3.write(list(indirect))
        # col3.write(f"\n\nALL associations: {len(list(all_associations))}")
        # col3.write(list(all_associations))

    else:
        st.write("Please select at least one phenotype term!")
    return None
