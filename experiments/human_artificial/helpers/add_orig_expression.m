%%% This script takes the simulation gene structure and adds the expression values to the transcripts for rquant evaluation

genes = '5000';
reads = '5000000';
chrms = '';
%stage = '2';

infile_matlab = sprintf('annotation/hg19_%ssubsample_%s_genes.nochr.mat', chrms, genes);
infile_fluxsim = sprintf('%s_genes_%s_reads/hg19_%ssubsample_%s_genes.gtf.pro', genes, reads, chrms, genes);
outfile = sprintf('%s_genes_%s_reads/genes.mat', genes, reads)

%%% load genes in matlab structure
load(infile_matlab);


%%% get expression information from the respective FluxSim pro file
data = importdata(infile_fluxsim);

for g = 1:length(genes),
    genes(g).expr_orig = [];
    for e = 1:length(genes(g).transcripts),
        e_idx = strmatch(genes(g).transcripts{e}, data.textdata(:, 2));
        if length(e_idx) == 1,
            genes(g).expr_orig(end + 1) = data.data(e_idx, 3);
        else
            genes(g).expr_orig(end + 1) = 0;
        end;
    end;
end;

save(outfile, 'genes');
%unix(['cd ' setsize '_reads/; ln -s ' strrep(anno_file, '.mat', '.expr.mat') ' genes.mat; cd ..']);

