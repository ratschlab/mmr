function genes = add_orig_expression(filename_matlab, filename_flux)

%%% This script takes the simulation gene structure and adds the expression values to the transcripts for rquant evaluation

%%% load genes in matlab structure
load(filename_matlab);

%%% get expression information from the respective FluxSim pro file
[names expression]= textread(filename_flux, '%*s%s%*s%*s%*s%d%*s%*s%*s%*s%*s%*s%*s');

for g = 1:length(genes),
    genes(g).expr_orig = [];
    for e = 1:length(genes(g).transcripts),
        e_idx = strmatch(genes(g).transcripts{e}, names);
        if length(e_idx) == 1,
            genes(g).expr_orig(end + 1) = expression(e_idx);
        else
            genes(g).expr_orig(end + 1) = 0;
        end;
    end;
end;

save(filename_matlab, 'genes');
