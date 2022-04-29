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
        # Mapping from term ID to name
        id_to_name = {id_: data.get("name") for id_, data in graph.nodes(data=True)}

        res = {key + " : " + id_to_name[key]: key for key in id_to_name.keys()}

        gene2hpo = {}
        # gene2hpo_assoc = {}
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
                gene2hpo[gene_name].append(hp_term)
                gene2hpo[gene_name] = list(
                    dict.fromkeys(gene2hpo[gene_name])
                )  # check for duplicate values and drop them
                gene2hpo[gene_name] = list(filter(None, gene2hpo[gene_name]))

                # gene2hpo_assoc[gene_name] = list(
                #    dict.fromkeys(gene2hpo_assoc[gene_name])
                # )  # check for duplicate values and drop them
                # gene2hpo_assoc[gene_name] = list(
                #    filter(None, gene2hpo_assoc[gene_name])
                # )

        gene2hpo_dict = gene_df(gene2hpo,graph)
        return gene2hpo, gene2hpo_dict, graph, id_to_name, res  # ,gene2hpo_assoc

    def get_graph(graph, hpos):
        final_graph = nx.MultiGraph()
        for term in hpos:
            sub_graph = nx.subgraph(graph, nx.predecessor(graph, term))
            final_graph = nx.compose(sub_graph, final_graph)
        final_graph = nx.DiGraph(final_graph)
        return final_graph

    @st.cache(allow_output_mutation=True)
    def gene_df(gene2hpo,graph):
        gene2hpo_dict = {}
        for gene_name in gene2hpo.keys():
            for hp_term in gene2hpo[gene_name]:
                if gene_name not in gene2hpo_dict:
                    gene2hpo_dict[gene_name] = list([])
                gene2hpo_dict[gene_name].append(hp_term)
                gene2hpo_dict[gene_name] = gene2hpo_dict[gene_name] + list(nx.predecessor(graph, hp_term))
            gene2hpo_dict[gene_name] = list(
                    dict.fromkeys(gene2hpo_dict[gene_name])
                )  # check for duplicate values and drop them

        gene2hpo_dict = {"Genes": list(gene2hpo_dict.keys()), "HPOs": list(gene2hpo_dict.values())}
        gene2hpo_dict = pd.DataFrame(gene2hpo_dict)
        return gene2hpo_dict

    @st.cache(allow_output_mutation=True)
    def gene_ranks(hpo_terms, gene_scores):
        term_len = len(hpo_terms)
        gene_scores["score"] = [
                            (len(list(set(hpo_terms) & set(i)))) / term_len
                            for i in gene_scores["HPOs"]
                        ]
        gene_scores = gene_scores.sort_values(by="score", ascending=False).reset_index(
                                drop=True
                            )

        gene_scores.drop("HPOs", axis=1, inplace=True)
        gene_scores.index = gene_scores.index + 1
        gene_scores = gene_scores.style.format(subset="score", precision=2).bar(subset="score", align="mid")

        return gene_scores

    data_load_state = st.text("Loading data...")
    gene2hpo, gene2hpo_dict, graph, id_to_name, res = load_data()
    data_load_state.text("Loaded HPO data!")
    data_load_state.text("")

    col1, col2, col3 = st.columns(3)

    # if col1.checkbox("Show loaded HPO metadata data"):
    #    col1.write(data.graph)
    terms = []
    #hpo_entry = col1.radio("Type of HPO input", ['List', 'Manual'])
    #if hpo_entry == 'Manual':
    #    options = col1.multiselect("Select HPO terms:", list(res.keys()))
    #    if col1.button('Submit'):
    #            terms = [res[term] for term in options]
    #else:
    #    hpo_text = col1.text_input('List of HPO terms seperated by comma/space/semi-colon', '')
    #    if hpo_text != '':
    #        terms = re.split(r'[,;\s]',hpo_text)

    options = col1.multiselect("Select HPO terms:", list(res.keys()))
    terms = [res[term] for term in options]
    hpo_text = col1.text_input('List of HPO terms seperated by comma/space/semi-colon', '')
    terms =  terms + re.split(r'[,;\s]',hpo_text)
    terms = list(filter(None, terms))

    gene = col3.selectbox("Select a Gene:", list(gene2hpo.keys()))
    col2.write("Hazel scores for Genes:")

    #if col1.button('Submit') and len(terms) != 0:
    if len(terms) != 0:
        term_len = len(terms)
        list_terms = {key + " : " + id_to_name[key] for key in terms}
        col1.write(f"Total HPO terms = {term_len}")
        col1.write(list(list_terms))

        gene_scores = gene_ranks(terms, gene2hpo_dict)
        #col2.dataframe(gene_scores)
        col2.write(gene_scores.to_html(), unsafe_allow_html=True)

        final_graph = get_graph(graph, gene2hpo[gene])

        gene2hpo_assoc = []
        for term in gene2hpo[gene]:
            gene2hpo_assoc = gene2hpo_assoc + list(nx.ancestors(final_graph.reverse(), term))
        gene2hpo_assoc = list(dict.fromkeys(gene2hpo_assoc))
        gene2hpo_assoc = list(set(gene2hpo_assoc) - set(gene2hpo[gene]))

        score = (len(list(set(terms) & set(gene2hpo[gene] + gene2hpo_assoc)))) / term_len
        col3.subheader(f"Hazel score = {score}")

        direct = {key + " : " + id_to_name[key] for key in list(set(terms) & set(gene2hpo[gene]))}
        indirect = {key + " : " + id_to_name[key] for key in list(set(terms) & set(gene2hpo_assoc))}
        #all_associations = {
        #    key + " : " + id_to_name[key]
        #    for key in list(set(terms) & set(gene2hpo[gene] + gene2hpo_assoc))
        #}

        col3.write(f"\n\nDirect associations: {len(list(direct))}")
        col3.write(list(direct))
        col3.write(f"\n\nIndirect associations: {len(list(indirect))}")
        col3.write(list(indirect))
        #col3.write(f"\n\nALL associations: {len(list(all_associations))}")
        #col3.write(list(all_associations))

    else:
        st.write("Please select at least one phenotype term!")
    return None
