def PlotEvidencesDistribution(table: pd.DataFrame, outfile: str = None):
    """
    This function takes a pandas dataframe with columns 'METABOLITE_id', 'reaction_EC' and an output file name as input,
    and groups the evidence codes by metabolite, counting the frequency of each evidence code per metabolite.
    It then plots a bar graph of the result and saves the figure in the output file (if outfile is not None)
    Args:
    table: a pandas dataframe with columns 'METABOLITE_id', 'reaction_EC'
    outfile: output file name, a string
    Returns:
    Plotly figure object
    """
    # Initialize empty dataframe
    nT = pd.DataFrame()

    # iterate over the rows of the input table
    for i, index in enumerate(table.index):
        # split the evidence codes by ;
        if isinstance(table.reaction_EC[index], str):
            evids = table.reaction_EC[index].split(";")
            for e in evids:
                if e in ecs_from_Uniprot:
                    nT.loc[i, "METABOLITE_id"] = table.METABOLITE_id[index]
                    nT.loc[i, "reaction_EC"] = f"{e}_HoGM"
                else:
                    nT.loc[i, "METABOLITE_id"] = table.METABOLITE_id[index]
                    nT.loc[i, "reaction_EC"] = e
        else:
            nT.loc[i, "METABOLITE_id"] = table.METABOLITE_id[index]
            nT.loc[i, "reaction_EC"] = "NoEC"
            i += 1

    # group and count the frequency of the evidence codes per metabolite
    group = (
        nT.groupby(["METABOLITE_id", "reaction_EC"])
        .size()
        .reset_index(name="Frequency")
    )
    group["Frequency"] = [1 for i in group.Frequency]

    # count total number of evidence codes
    group["EC_count"] = group.reaction_EC.map(
        group.groupby("reaction_EC")["Frequency"].count().to_dict()
    )

    # sort the dataframe by EC_count
    group.sort_values("EC_count", ascending=False, inplace=True, ignore_index=True)

    # plot the graph
    fig = px.bar(
        data_frame=group,
        y="Frequency",
        x="reaction_EC",
        color="METABOLITE_id",
        orientation="v",
        pattern_shape="METABOLITE_id",
        color_discrete_sequence=px.colors.qualitative.Pastel,
        pattern_shape_sequence=["", "/", "\\", "x", "", "-", "|", "+", ".", ""],
    ).update_xaxes(categoryorder="total descending")

    # update graph layout
    fig.update_layout(
        title_text="Evidences Per Result",
        height=600,
        width=1200,
        plot_bgcolor="rgba(0, 0, 0, 0)",
        font=dict(size=12),
    )
    if outfile:
        fig.write_image(outfile, height=600, width=1200, validate=True, scale=2)
    else:
        pass
    return fig


def PlotMetaboliteInfo(singleMetTable: pd.DataFrame, outfile: str = None):
    """
    This function takes a pandas dataframe  with information of a single metabolite, and an output file name as input.
    It creates a subplot with two rows and two columns with the following information:
    1. Table of metabolite information from the input dataframe
    2. Plot of isotopic distribution (bar chart) of the metabolite
    3. Image of the molecular structure of the metabolite
    The function saves the figure in the output file (if outfile is not None)

    Args:
    singleMetTable: a pandas dataframe with information of a single metabolite
    outfile: output file name, a string. Default None
    Returns:
    Plotly figure object
    """

    df = singleMetTable.T.reset_index()
    mol = Chem.MolFromSmiles(singleMetTable.METABOLITE_smiles.iloc[0])
    d2d = rdMolDraw2D.MolDraw2DCairo(400, 200)
    d2d.ClearDrawing()
    dos = d2d.drawOptions()
    dos.bondLineWidth = 3
    dos.centreMoleculesBeforeDrawing = True
    dos.useBWAtomPalette()
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    img = Image.open(io.BytesIO(d2d.GetDrawingText()))

    mz = {
        float(value.split(":")[0]): float(value.split(":")[1])
        for value in singleMetTable.MZ_distribution.iloc[0].split(";")
    }

    F = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "table"}, {"type": "bar"}]],
        subplot_titles=("Result information", "Isotopic distribution"),
    )

    F.add_trace(
        go.Bar(name="Metabolite", x=tuple(mz.keys()), y=tuple(mz.values())),
        row=1,
        col=2,
    )

    F.add_trace(
        go.Table(
            header=dict(values=["DATA", "VALUE"], font=dict(size=12), align="left"),
            cells=dict(values=[df[k].tolist() for k in df.columns], align="left"),
        ),
        row=1,
        col=1,
    )

    F.add_layout_image(
        dict(
            source=img,
            x=0.95,
            y=0.5,
        ),
        row=1,
        col=2,
    )

    F.update_layout_images(
        dict(
            xref="paper",
            yref="paper",
            sizex=0.3,
            sizey=0.35,
            xanchor="right",
            yanchor="bottom",
        )
    )

    F.update_traces(width=0.05, row=1, col=2)
    F.update_yaxes(range=[0, 100], row=1, col=2)
    F.update_layout(
        width=1200,
        height=600,
        plot_bgcolor="rgba(0, 0, 0, 0)",
        font=dict(size=12),
        xaxis_title="Atomic Mass (u)",
        yaxis_title="Relative abundance (%)",
    )
    if outfile:
        F.write_image(
            outfile, height=400 * rows, width=400 * cols, validate=True, scale=2
        )
    else:
        pass
    return F


def PlotClusterMap(
    singleMetabolite: pd.DataFrame,
    outfile: str = None,
    selection: str = "all",
    n_samples: int = 75,
):
    """
    This function plots a heatmap of pairwise alignments of sequences of the metabolite, and clusters the sequences based on
    their identity percentage.

    Parameters:
    singleMetabolite (pd.DataFrame): A dataframe containing information of the metabolite of interest.
    outfile (str, optional): A filepath to save the plot.
    selection (list, optional): A list of evidence codes to filter the sequences.
    n_samples (int, optional): The number of sequences to sample.

    Returns:
    None

    """
    if selection == "all":
        sample = singleMetabolite.copy()
    else:
        sample = singleMetabolite[singleMetabolite.EVIDENCE.isin(selection)]

    print(
        f"{len(sample.index)} unique sequences found for {len(sample.ORGANISM_name.unique())} species"
    )
    if len(sample.index) >= n_samples:
        print(
            f"Sampling ~{n_samples}sequences with higher existence and species diversity"
        )
        max_quality = sample[sample.EXISTENCE == "Evidence at protein level"]
        max_quality = max_quality.groupby(["EVIDENCE", "ORGANISM_name"]).apply(
            lambda g: g.sample(int(1 * len(g) / len(sample) + 1))
        )
        if len(max_quality) < n_samples:
            sample = sample[sample.EXISTENCE != "Evidence at protein level"]
            sample = sample.groupby(["EVIDENCE"]).apply(
                lambda g: g.sample(
                    int((n_samples - len(max_quality)) * len(g) / len(sample) + 1)
                )
            )
            sample = pd.concat([sample, max_quality], ignore_index=True)
            sample.sort_values(by="EXISTENCE", ignore_index=True, inplace=True)
            sample.drop_duplicates(ignore_index=True, inplace=True)
        else:
            sample = max_quality.sample(n_samples, ignore_index=True)
    else:
        pass

    hmap = pd.DataFrame()
    for index in sample.index:
        for jndex in sample.index:
            a = sample.SEQUENCE[index]
            b = sample.SEQUENCE[jndex]
            alignment = pairwise2.align.globalxx(a, b, score_only=True)
            identity = (alignment * 100) / len(b)
            hmap.loc[
                f'{sample.ENTRY_id[index]} : {sample.ORGANISM_name[index].split("(")[0]}',
                sample.ENTRY_id[jndex],
            ] = identity

    lut_phylum = dict(zip(set(sample.SUPERKINGDOM), sns.color_palette("Set1", 10)))
    row_colors = sample.SUPERKINGDOM.map(lut_phylum)

    lut_ev = dict(zip(set(sample.EVIDENCE), sns.color_palette("Set3", 15)))
    row_colors2 = sample.EVIDENCE.map(lut_ev)

    lut_ex = dict(zip(set(sample.EXISTENCE), sns.color_palette("Set2", 15)))
    col_colors = sample.EXISTENCE.map(lut_ex)

    ax = sns.clustermap(
        hmap.reset_index(drop=True),
        figsize=(17, 15),
        xticklabels=1,
        yticklabels=hmap.index,
        cmap="RdBu_r",
        method="complete",
        vmin=0,
        vmax=100,
        tree_kws=dict(linewidths=1.5),
        row_colors=[row_colors, row_colors2],
        col_colors=[col_colors],
    )

    x0, _y0, _w, _h = ax.cbar_pos
    ax.ax_cbar.set_position([x0 + 0.05, 0.82, _w / 2, _h / 1.3])
    ax.ax_cbar.set_ylabel("% Identity", fontsize=16, fontweight="bold")
    ax.ax_cbar.tick_params(axis="y", length=3, width=1, labelsize=12)

    for label in sample.SUPERKINGDOM.unique():
        ax.ax_row_dendrogram.bar(
            0, 0, color=lut_phylum[label], label=label, linewidth=0
        )
    l1 = ax.ax_row_dendrogram.legend(
        title="Taxonomy",
        loc="center",
        ncol=1,
        bbox_to_anchor=(0.05, 0.03),
        bbox_transform=plt.gcf().transFigure,
        frameon=False,
    )

    for label in sample.EXISTENCE.unique():
        ax.ax_col_dendrogram.bar(0, 0, color=lut_ex[label], label=label, linewidth=0)
    l2 = ax.ax_col_dendrogram.legend(
        title="Existence",
        loc="center",
        ncol=1,
        bbox_to_anchor=(0.87, 0.82),
        bbox_transform=plt.gcf().transFigure,
        frameon=False,
    )

    xx = []
    for label in sample.EVIDENCE.unique():
        x = ax.ax_row_dendrogram.bar(
            0, 0, color=lut_ev[label], label=label, linewidth=0
        )
        xx.append(x)

    l3 = plt.legend(
        xx,
        sample.EVIDENCE.unique(),
        loc="center",
        title="Evidence",
        bbox_to_anchor=(0.15, 0.03),
        bbox_transform=plt.gcf().transFigure,
        frameon=False,
    )

    if isinstance(outfile, str):
        ax.savefig(outfile, dpi=300, format="svg", bbox_inches="tight")
    return sample, ax


def PlotSequenceLogo(
    table: pd.DataFrame,
    prefix: str = "MSA",
    selection: str = "all",
    outfile: str = None,
) -> None:
    """
    Generate sequence logo for a multiple sequence alignment of domain sequences.

    Parameters:
    table (pd.DataFrame): A DataFrame containing information of domains' location and sequences.
    prefix (str, optional): The prefix of the multiple sequence alignment file. Defaults to 'MSA'.
    selection (str, optional): The option to select domains to be included in the multiple sequence alignment. Defaults to 'all'.
    outfile (str, optional): The output file path of the sequence logo plot. Defaults to None.

    Returns:
    None
    """

    GenerateMSA(table=table, prefix=prefix, selection=selection)

    domsTable = (
        table.dropna(subset=["DOMAIN_location"])[
            ["ENTRY_id", "DOMAIN_location", "LENGHT"]
        ]
        .sort_values(by="LENGHT", ascending=True)
        .reset_index(drop=True)
    )

    domain_pattern = re.compile("\d+\.\.\d+")
    names_pattern = re.compile('note="([^"]*)"')

    domains_dict = {}
    entries = {}
    for index in domsTable.index:
        namesDomains = names_pattern.findall(domsTable.DOMAIN_location[index])
        for didx, name in enumerate(namesDomains):
            domains_dict[name] = list(
                map(
                    int,
                    domain_pattern.findall(domsTable.DOMAIN_location[index])[
                        didx
                    ].split(".."),
                )
            )
            entries[name] = domsTable.ENTRY_id[index]

    alignment = AlignIO.read(f"data/MSA/{prefix}.aln", "fasta")

    msa_dict = {}
    for name, entry in entries.items():
        for aln in alignment:
            if entry in aln.id:
                new_numeration = list(
                    {i + 1: L for i, L in enumerate(aln.seq) if L != "-"}.items()
                )
                start, end = domains_dict[name]
                msa_dict[name] = [
                    new_numeration[start - 1][0],
                    new_numeration[end - 1][0],
                ]

    chunks = {}
    for seq in alignment:
        sequence_chunks = Chunk_list(str(seq.seq), 100)
        for index, seq_fragment in enumerate(sequence_chunks):
            chunks.setdefault(index, []).append(seq_fragment)

    fig, axes = plt.subplots(
        nrows=len(chunks),
        ncols=1,
        figsize=(20, 2 * len(chunks)),
        sharex=False,
        sharey=False,
    )
    xtick_start = 1

    yticks = np.arange(0, len(alignment) + 1, 5)

    for index, seqs in chunks.items():
        ax = axes.ravel()[index]
        counts_mat = logomaker.alignment_to_matrix(
            sequences=seqs,
            to_type="counts",
        )

        logo = logomaker.Logo(
            counts_mat,
            shade_below=0.5,
            fade_below=0.5,
            ax=ax,
            color_scheme="black",
            fade_probabilities=True,
            stack_order="big_on_top",
        )

        logo.style_spines(visible=False)
        logo.style_xticks(rotation=90, fmt="%d", anchor=0, fontsize=10)
        logo.style_spines(spines=["left", "bottom"], visible=True)
        ax.set_ylabel("Bits", labelpad=1, size=16)
        # ax.set_yticklabels(ax.get_yticks().astype(int),fontsize=10)

        xticklabels = ax.get_xticks() + xtick_start

        ax.set_xticklabels(xticklabels)
        ax.tick_params(width=1.5, size=2)

        colorsDomains = {
            k: color
            for k, color in zip(
                msa_dict.keys(),
                sns.color_palette("Set1", len(msa_dict.keys())).as_hex(),
            )
        }

        for dname, msa_pointers in msa_dict.items():
            s = msa_pointers[0]
            e = msa_pointers[1]
            if s in xticklabels:
                plot_index = list(xticklabels).index(s)
                ax.axvline(
                    x=ax.get_xticks()[plot_index],
                    ymin=0,
                    ymax=len(alignment) + 1,
                    color=colorsDomains[dname],
                    linewidth=6,
                    linestyle="--",
                )
                ax.text(
                    x=ax.get_xticks()[plot_index] + 0.5,
                    y=len(alignment),
                    s=f"{dname}->",
                    fontsize=12,
                    horizontalalignment="left",
                    fontweight="bold",
                    color=colorsDomains[dname],
                )
            if e in xticklabels:
                plot_index = list(xticklabels).index(e)
                ax.axvline(
                    x=ax.get_xticks()[plot_index],
                    ymin=0,
                    ymax=len(alignment) + 1,
                    color=colorsDomains[dname],
                    linewidth=6,
                    linestyle="--",
                )
                ax.text(
                    x=ax.get_xticks()[plot_index] - 0.5,
                    y=len(alignment),
                    s=f"<-{dname}",
                    fontsize=12,
                    horizontalalignment="right",
                    fontweight="bold",
                    color=colorsDomains[dname],
                )

        xtick_start = xtick_start + ax.get_xticks()[-1] + 1

    for ax in axes:
        ax.set_ylim(0, len(alignment) + 1)
        ax.set_yticks(yticks, yticks)

    fig.subplots_adjust(
        bottom=0.0, top=1.0, left=0.1, right=0.8, wspace=0.2, hspace=0.4
    )

    if isinstance(outfile, str):
        fig.savefig(outfile, dpi=300, format="svg", bbox_inches="tight")

    else:
        pass

    return fig


def PlotGO(table: pd.DataFrame, outfile: str = None, selection: str = "all"):
    """
    PlotGO function is used for visualizing the Gene ontology (GO) of the given table with bar graph.
    :param table: pandas dataframe which contains the information of the GO terms
    :param outfile: a string which represents the path to save the output file, if None the plot will only be shown.
    :param selection: a string which represents the type of evidence, default is 'all'
    """
    if selection == "all":
        table = table
    else:
        table = table[table.EVIDENCE.isin(selection)]

    F = make_subplots(
        rows=1,
        cols=6,
        shared_xaxes=False,
        shared_yaxes=False,
        horizontal_spacing=0,
        subplot_titles=["GO_process", "", "", "", "", "GO_function"],
    )
    go_table = pd.DataFrame()
    for index in table.index:
        row = table.loc[index]
        if isinstance(row.GO_process, str):
            GO_terms = row.GO_process.split(";")
            for go in GO_terms:
                row.ORGANISM_name = row.ORGANISM_name.split("(")[0]
                row.GO_process = go.strip().split(",")[-1]
                go_table = go_table.append(row, ignore_index=True)
        else:
            row.GO_process = "Not available"
            go_table = go_table.append(row, ignore_index=True)

    go_process = (
        go_table.groupby(["GO_process", "ORGANISM_name"])["SEQUENCE"]
        .count()
        .reset_index(["GO_process", "ORGANISM_name"])
    )

    go_process.sort_values(by="ORGANISM_name", inplace=True, ignore_index=True)

    fig1 = px.bar(
        data_frame=go_process,
        x="SEQUENCE",
        y="GO_process",
        color="ORGANISM_name",
        orientation="h",
        color_discrete_sequence=px.colors.qualitative.Pastel,
    )

    go_table = pd.DataFrame()
    for index in table.index:
        row = table.loc[index]
        if isinstance(row.GO_function, str):
            GO_terms = row.GO_function.split(";")
            for go in GO_terms:
                row.ORGANISM_name = row.ORGANISM_name.split("(")[0]
                row.GO_function = go.strip().split(",")[-1]
                go_table = go_table.append(row, ignore_index=True)
        else:
            row.GO_function = "Not available"
            go_table = go_table.append(row, ignore_index=True)

    go_function = (
        go_table.groupby(["GO_function", "ORGANISM_name"], dropna=False)["SEQUENCE"]
        .count()
        .reset_index(["GO_function", "ORGANISM_name"])
    )

    go_function.sort_values(by="ORGANISM_name", inplace=True, ignore_index=True)

    fig2 = px.bar(
        data_frame=go_function,
        x="SEQUENCE",
        y="GO_function",
        color="ORGANISM_name",
        orientation="h",
        color_discrete_sequence=px.colors.qualitative.Pastel,
    )

    for d in fig1.data:
        d.showlegend = False
        F.add_trace(d, row=1, col=1)

    for d in fig2.data:
        d.showlegend = True
        F.add_trace(d, row=1, col=6)

    F.update_xaxes(
        title_text="Frequency",
        showgrid=True,
        row=1,
        col=1,
    )
    F.update_xaxes(title_text="Frequency", showgrid=True, row=1, col=6)

    F.update_layout(
        barmode="stack",
        height=800,
        width=1200,
        title_text=f"GO Associations (Organism names by UniProt)",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        font=dict(size=11),
    )
    if outfile:
        F.write_image(outfile, height=800, width=1200, validate=True, scale=2)
    else:
        pass
    return F
