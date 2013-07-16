function eval_correlation_bowtie(readlen)
%function eval_correlation_bowtie(readlen)

    %%% import corralation function mycorr
    addpath ~/git/projects/2011/rquant/evaluation/

    readnum = '500000';
    gene_num = '5000';

    noise_levels = {''};%, 'noise0.01', 'noise0.02', 'noise0.03'};
    %experiments = {'unfiltered', 'best', 'mmr0'};
    experiments = {'unfiltered', 'best', 'mmr0', 'mmr1', 'mmr2', 'mmr3', 'mmr4', 'mmr5'};
    %experiments = {'unfiltered', 'best', 'mmr0', 'mmr1', 'mmr2'};

    experiment = [gene_num '_genes_' readnum '_reads'];

    %chrms = 'chr2_chr3_ch4_';
    chrms = '';

    %%% do not load artifical data, but instead assume, that the gene structure carries the field 'expr_orig'

    pearson = [];
    spearman = [];
    labels = {};

    for n_idx = 1:length(noise_levels),

        noise = noise_levels{n_idx};
        labels{end + 1} = [' Bowtie2 ' noise];

        for e_idx = 1:length(experiments),

            which_set = experiments{e_idx};

            fprintf(1, 'Evaluating %s %s\n====================\n\n', which_set, noise);

            %%% load art data into genes
            load(sprintf('%s/genes.mat', experiment));
            if ~isfield(genes, 'expr_orig'),
                genes = add_orig_expression(sprintf('%s/genes.mat', experiment), sprintf('%s_genes_%s_reads/hg19_%ssubsample_%s_genes.gtf.pro', gene_num, readnum, chrms, gene_num));
            end;
            if ~isfield(genes, 'transcript_length'),
                genes = add_transcript_length(sprintf('%s/genes.mat', experiment));
            end;

            art_orig = genes;

            %%% extract all transcripts and their IDs from the gene structure
            art.IDs = [art_orig(:).transcripts]';
            art.expr = [art_orig(:).expr_orig]';
            art.lengths = [art_orig(:).transcript_length]';

            %%% sort by transcript ID
            [art.IDs, s_idx] = sort(art.IDs);
            art.expr = art.expr(s_idx);
            art.lengths = art.lengths(s_idx);

            %%% load bowtie results
            f_tag = ['.' which_set];
            if strcmp(f_tag, '.unfiltered')
                f_tag = '';
            end;
            n_tag = ['.' noise];
            if strcmp(n_tag, '.'),
                n_tag = '';
            end;
            [a b] = textread(sprintf('%s/bowtie/hg19_%ssubsample_%s_genes%s/alignment%s.counts', experiment, chrms, gene_num, n_tag, f_tag), '%f%s');
            if ~strcmp(which_set, 'unfiltered'),
                a = a(2:end);
                b = b(2:end);
            end;
            pred.textdata = b;
            pred.data = a;

            %%% remove NaNs
            s_idx = find(isnan(pred.data));
            pred.data(s_idx) = [];
            pred.textdata(s_idx) = [];

            %%% sort by transcript ID
            [tmp s_idx] = sort(pred.textdata);
            pred.textdata = pred.textdata(s_idx);
            pred.data = pred.data(s_idx);

            %%% intersect artificial data and prediction
            [tmp s1 s2] = intersect(art.IDs, pred.textdata);
            art.IDs = art.IDs(s1);
            art.expr = art.expr(s1);
            art.lengths = art.lengths(s1);
            pred.textdata = pred.textdata(s2);
            pred.data = pred.data(s2);

            r_p = mycorr(art.expr, pred.data .* readlen ./ art.lengths, 'Pearson');
            r_s = mycorr(art.expr, pred.data .* readlen ./ art.lengths, 'Spearman');
            pearson(n_idx, e_idx) = r_p;
            spearman(n_idx, e_idx) = r_s;
            fprintf(1, 'Bowtie (Pearson / Spearman):\n\ttranscripts: %i (%i)\n\tr: %f / %f\n\n', size(art.expr, 1), size([art_orig(:).transcripts]', 1), r_p, r_s);
        end;
    end;

    clf;
    colormap('summer');

    subplot(1, 2, 1);
    bar(pearson);
    title('Bowtie (Pearson)');
    set(gca, 'XTick', 1:length(noise_levels));
    noise_labels = {'native', '+1%', '+2%', '+3%'};
    set(gca, 'XTickLabel', noise_labels(1:length(noise_levels)));
    xlabel('Error Rate');
    ylabel('Transcript Correlation (Pearson)');
    ylim([0.0 1.0]);
    %legend(experiments, 'Location', 'SouthEast');

    subplot(1, 2, 2);
    bar(spearman);
    title('Bowtie (Spearman)');
    set(gca, 'XTick', 1:length(noise_levels));
    noise_labels = {'native', '+1%', '+2%', '+3%'};
    set(gca, 'XTickLabel', noise_labels(1:length(noise_levels)));
    xlabel('Error Rate');
    ylabel('Transcript Correlation (Spearman)');
    ylim([0.0 1.0]);
    legend(experiments, 'Location', 'NorthEast');

    print('-dpdf', '-S1200,600', sprintf('%s/transcript_correlation_bowtie.pdf', experiment));

    return

