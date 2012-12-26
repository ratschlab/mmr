
readnum = '15000000';
gene_num = '5000';

%noise_levels = {'noise0.01', 'noise0.02', 'noise0.03', ''};
noise_levels = {''};
experiments = {'unfiltered', 'best', 'mmr0'};

experiment = [gene_num '_genes_' readnum '_reads'];

%chrms = 'chr2_chr3_ch4_';
chrms = '';

%%%% load artificial data
%%%% art.data -> 1 length, 2 expr fraction, 3 expr number, 4 libr fraction, 5 libr number, 6 seq fraction, 7 seq number, 8 cov fraction, 9 chi square, 10 coeff of var
%%%% art.textdata -> 1 locus, 2 trans ID, 3 coding (CDS|NC)
%sprintf('%s/hg19_%ssubsample_%s_genes.gtf.pro', experiment, chrms, gene_num)
%art_orig = importdata(sprintf('%s/hg19_%ssubsample_%s_genes.gtf.pro', experiment, chrms, gene_num));

%%% do not load artifical data, but instead assume, that the gene structure carries the field 'expr_orig'

pearson_cuff = [];
pearson_rquant = [];
spearman_cuff = [];
spearman_rquant = [];
labels_rquant = {};
labels_cuff = {};

for n_idx = 1:length(noise_levels),

    noise = noise_levels{n_idx};
    labels_cuff{end + 1} = [' Cufflinks ' noise];
    labels_rquant{end + 1} = [' rQuant ' noise];

    for e_idx = 1:length(experiments),

        which_set = experiments{e_idx};

        fprintf(1, 'Evaluating %s %s\n====================\n\n', which_set, noise);

        %%% load art data into genes
        load(sprintf('%s/genes.mat', experiment));
        if ~isfield(genes, 'expr_orig'),
            error(sprintf('Gene structure %s/genes.mat lacks expression information\n', experiment));
        end;
        if ~isfield(genes, 'multi_frac_tophat'),
            genes = add_multimap_info_tophat(sprintf('%s/genes.mat', experiment), sprintf('%s/tophat/hg19_subsample_%s_genes.gtf.fastq.gz/accepted_hits', experiment, gene_num));
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

        %%% load cufflinks results
        f_tag = ['.' which_set];
        if strcmp(f_tag, '.unfiltered')
            f_tag = '';
        end;
        n_tag = [noise '.'];
        if strcmp(n_tag, '.'),
            n_tag = '';
        end;
        pred = importdata(sprintf('%s/tophat/hg19_%ssubsample_%s_genes.gtf.%sfastq.gz/cufflinks/%s/transcripts.counts', experiment, chrms, gene_num, n_tag, which_set));
        pred.textdata = pred.textdata(2:end, :); %% textdata -> 1 genes, 2 transcripts; data -> 1 FPKM, 2 cov

        %%% remove NaNs
        s_idx = find(isnan(pred.data(:, 1)));
        pred.data(s_idx, :) = [];
        pred.textdata(s_idx, :) = [];

        %%% sort by transcript ID
        [tmp s_idx] = sort(pred.textdata(:, 2));
        pred.textdata = pred.textdata(s_idx, :);
        pred.data = pred.data(s_idx, :);

        %%% intersect artificial data and prediction
        [tmp s1 s2] = intersect(art.IDs, pred.textdata(:, 2));
        art.IDs = art.IDs(s1);
        art.expr = art.expr(s1);
        art.lengths = art.lengths(s1);
        pred.textdata = pred.textdata(s2, :);
        pred.data = pred.data(s2, :);

        %%% TODO Cufflinks uses rpkms
        [r_p p_p] = corr(art.expr, pred.data(:, 1) .* art.lengths ./ 1000, 'type', 'Pearson');
        [r_s p_s] = corr(art.expr, pred.data(:, 1) .* art.lengths ./ 1000, 'type', 'Spearman');
        pearson_cuff(n_idx, e_idx) = r_p;
        spearman_cuff(n_idx, e_idx) = r_s;
        %fprintf(1, 'Cufflinks (Pearson / Spearman):\n\ttranscripts: %i (%i)\n\tr: %f / %f\n\tp:%f / %f\n\n', size(art.data, 1), size(art_orig.data, 1), r_p, r_s, p_p, p_s);
        fprintf(1, 'Cufflinks (Pearson / Spearman):\n\ttranscripts: %i (%i)\n\tr: %f / %f\n\n', size(art.expr, 1), size([art_orig(:).transcripts]', 1), r_p, r_s);

        %for perc = [0.05 0.1 0.15 0.2 0.25],
        for perc = [0.1 1.0],

            %%% extract all transcripts and their IDs from the gene structure
            art.IDs = [art_orig(:).transcripts]';
            art.expr = [art_orig(:).expr_orig]';
            art.multi_frac = [art_orig(:).multi_frac_tophat]';

            %%% filter by mm percentile
            [art.multi_frac s_idx] = sort(art.multi_frac, 1, 'descend');
            art.expr = art.expr(s_idx);
            art.IDs = art.IDs(s_idx);
                
            p = floor(length(art.multi_frac) * perc); 
            art.expr = art.expr(1:p);
            art.IDs = art.IDs(1:p);

            %%% sort by transcript ID
            [art.IDs, s_idx] = sort(art.IDs);
            art.expr = art.expr(s_idx);

            %%% load rQuant data
            d_tag = which_set;
            n_tag = [noise '.'];
            if strcmp(n_tag, '.'),
                n_tag = '';
            else
                d_tag = [d_tag '.' noise];
            end;
            load(sprintf('%s/tophat/hg19_%ssubsample_%s_genes.gtf.%sfastq.gz/rquant/%s/accepted_hits%s_rquant.mat', experiment, chrms, gene_num, n_tag, d_tag, f_tag));

            %%% extract transcript info
            pred = struct();
            pred.IDs = [genes(:).transcripts]';
            pred.expr = [genes(:).transcript_weights]';
            %%% get rid of NaNs
            s_idx = find(isnan(pred.expr));
            pred.expr(s_idx) = [];
            pred.IDs(s_idx, :) = [];
            %%% sort
            [pred.IDs s_idx] = sort(pred.IDs);
            pred.expr = pred.expr(s_idx);

            %%% intersect artificial data and prediction
            [tmp s1 s2] = intersect(art.IDs, pred.IDs);
            art.IDs = art.IDs(s1);
            art.expr = art.expr(s1);
            pred.IDs = pred.IDs(s2);
            pred.expr = pred.expr(s2);

            %keyboard;
            %figure; plot(log10(art.data(:, 3)), log10(pred.data), 'o');
            %title(which_set);

            [r_p p_p] = corr(art.expr, pred.expr, 'type', 'Pearson');
            [r_s p_s] = corr(art.expr, pred.expr, 'type', 'Spearman');
            pearson_rquant(n_idx, e_idx) = r_p;
            spearman_rquant(n_idx, e_idx) = r_s;
            %fprintf(1, 'rQuant (Pearson / Spearman):\n\ttranscripts: %i (%i)\n\tr: %f / %f\n\tp:%f / %f\n\n', size(art.data, 1), size(art_orig.data, 1), r_p, r_s, p_p, p_s);
            fprintf(1, 'rQuant (Pearson / Spearman) - percentile %f:\n\ttranscripts: %i (%i)\n\tr: %f / %f\n\n', perc, size(art.expr, 1), size([art_orig(:).transcripts]', 1), r_p, r_s);
        end;
    end;
end;

clf;
plot(pearson_rquant', '-o');
title('Per Transcript Correlation (Pearson) - rQuant');
legend(labels_rquant, 'Location', 'NorthWest');
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', experiments);
xlim([0 4]);
print('-dpdf', sprintf('%s/transcript_correlation_rquant_pearson.pdf', experiment));

clf;
plot(pearson_cuff', '-o');
title('Per Transcript Correlation (Pearson) - Cufflinks');
legend(labels_cuff, 'Location', 'NorthWest');
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', experiments);
xlim([0 4]);
print('-dpdf', sprintf('%s/transcript_correlation_cufflinks_pearson.pdf', experiment));

clf;
plot(spearman_rquant', '-o');
title('Per Transcript Correlation (Pearson) - rQuant');
legend(labels_rquant, 'Location', 'NorthWest');
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', experiments);
xlim([0 4]);
print('-dpdf', sprintf('%s/transcript_correlation_rquant_spearman.pdf', experiment));

clf;
plot(spearman_cuff', '-o');
title('Per Transcript Correlation (Spearman) - Cufflinks');
legend(labels_cuff, 'Location', 'NorthWest');
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', experiments);
xlim([0 4]);
print('-dpdf', sprintf('%s/transcript_correlation_cufflinks_spearman.pdf', experiment));
