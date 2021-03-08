from dna_features_viewer import GraphicFeature, GraphicRecord


def plot_oligo_design(design, ax=None, figure_width=15, plot_lines=1, colors=None):
    
    if colors is None:
        colors = {'adaptor' : '#beaed4',
                  'oligo' : '#fdc086',
                  'annealing' : '#ffff99',
                  'target' : '#7fc97f',}

    gfeatures = []

    gfeatures.append( GraphicFeature(start=design.seq.index(design._tgt_designer.adt_lft), 
                                     end=design.seq.index(design._tgt_designer.adt_lft)+len(design._tgt_designer.adt_lft),
                                     label='adaptor',
                                     strand=+1, color=colors['adaptor']) ) #adaptor_left
    for i in range(len(design.annealings)):
        #olg_strandness = +1 if i%2==0 else -1
        gfeatures.append( GraphicFeature(start=design.seq.index(design.oligos_on_pos_strand[i]), 
                                        end=design.seq.index(design.oligos_on_pos_strand[i])+len(design.oligos_on_pos_strand[i]),
                                        label='oligo',
                                        strand=+1, color=colors['oligo']) ) # oligo
        gfeatures.append( GraphicFeature(start=design.seq.index(design.annealings[i]), 
                                        end=design.seq.index(design.annealings[i])+len(design.annealings[i]),
                                        label='annealing',
                                        strand=+1, color=colors['annealing']) ) # annealing
    i+=1
    gfeatures.append( GraphicFeature(start=design.seq.index(design.oligos_on_pos_strand[i]), 
                                    end=design.seq.index(design.oligos_on_pos_strand[i])+len(design.oligos_on_pos_strand[i]),
                                    label='oligo',
                                    strand=+1, color=colors['oligo']) ) # the last oligo
    gfeatures.append( GraphicFeature(start=design.seq.index(design._tgt_designer.adt_rgt), 
                                     end=design.seq.index(design._tgt_designer.adt_rgt)+len(design._tgt_designer.adt_rgt),
                                     label='adaptor',
                                     strand=+1, color=colors['adaptor']) ) #adaptor_right

    gfeatures.append( GraphicFeature(start=design.tgt_offset, 
                                     end=design.tgt_offset+len(design._tgt_designer.tgt_seq),
                                     label='target',
                                     strand=+1, color=colors['target']) ) # target

    record = GraphicRecord(sequence=design.seq, features=gfeatures)

    if plot_lines>1:
        record.plot_on_multiple_lines(plot_lines, plot_sequence=True, figure_width=figure_width)
    else:
        record.plot(ax=ax, plot_sequence=True, figure_width=figure_width)

