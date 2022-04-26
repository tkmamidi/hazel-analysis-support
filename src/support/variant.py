import streamlit as st
import pandas as pd
import obonet
import requests
import networkx as nx

# st.sidebar.subheader("OMIM diseases and associated HPO terms!")


def variant():
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
                # gene2hpo[gene_name] = gene2hpo[gene_name] + list(nx.neighbors(graph, hp_term))
                # gene2hpo_assoc[gene_name] = gene2hpo[gene_name] + list(nx.ancestors(graph.reverse(), hp_term))
                # gene2hpo[gene_name] = gene2hpo[gene_name] + list(graph.predecessors(hp_term)) + list(graph.successors(hp_term))
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

        return gene2hpo, graph, id_to_name, res  # ,gene2hpo_assoc

    def get_graph(graph, hpos):
        final_graph = nx.MultiGraph()
        for term in hpos:
            sub_graph = nx.subgraph(graph, nx.predecessor(graph, term))
            final_graph = nx.compose(sub_graph, final_graph)
        final_graph = nx.DiGraph(final_graph)
        return final_graph

    data_load_state = st.text("Loading data...")
    gene2hpo, graph, id_to_name, res = load_data()
    data_load_state.text("Loaded HPO data!")
    data_load_state.text("")

    col1, col2 = st.columns(2)

    # if col1.checkbox("Show loaded HPO metadata data"):
    #    col1.write(data.graph)

    options = col1.multiselect("Select HPO terms:", list(res.keys()))
    gene = col2.selectbox("Select a Gene:", list(gene2hpo.keys()))
    terms = [res[term] for term in options]

    if st.button("submit") and len(terms) != 0:
        term_len = len(terms)
        list_terms = {key + " : " + id_to_name[key] for key in terms}
        col1.write(f"Total HPO terms = {term_len}")
        col1.write(list(list_terms))

        final_graph = get_graph(graph, gene2hpo[gene])

        gene2hpo_assoc = []
        for term in gene2hpo[gene]:
            gene2hpo_assoc = gene2hpo_assoc + list(nx.ancestors(final_graph.reverse(), term))
        gene2hpo_assoc = list(dict.fromkeys(gene2hpo_assoc))
        gene2hpo_assoc = list(set(gene2hpo_assoc) - set(gene2hpo[gene]))

        score = (len(list(set(terms) & set(gene2hpo[gene] + gene2hpo_assoc)))) / term_len
        col2.subheader(f"Hazel score = {score}")

        direct = {key + " : " + id_to_name[key] for key in list(set(terms) & set(gene2hpo[gene]))}
        indirect = {key + " : " + id_to_name[key] for key in list(set(terms) & set(gene2hpo_assoc))}
        all_associations = {
            key + " : " + id_to_name[key]
            for key in list(set(terms) & set(gene2hpo[gene] + gene2hpo_assoc))
        }

        col2.write(f"\n\nDirect associations: {len(list(direct))}")
        col2.write(list(direct))
        col2.write(f"\n\nIndirect associations: {len(list(indirect))}")
        col2.write(list(indirect))
        col2.write(f"\n\nALL associations: {len(list(all_associations))}")
        col2.write(list(all_associations))

    return None
