function genes = add_multimap_info(filename, bam_name_stub)

    addpath '~/git/tools/ngs';

    load(filename);

    bam_mm = [bam_name_stub, '.multimappers.sorted.bam'];
    bam_unf = [bam_name_stub, '.sorted.bam'];

    for i = 1:length(genes),
        if mod(i, 10) == 0,
            fprintf(1, '.');
            if mod(i, 100) == 0,
                fprintf(1, '%i genes \n', i);
            end;
        end;
        gene = genes(i);
        mask1 = get_reads(bam_mm, gene.chr, gene.start, gene.stop, '0', 0, 1000, 50000, 2, 10);
        if isempty(mask1),
            reads1 = zeros(0,gene.stop-gene.start+1);
        else
            reads1 = sparse(mask1(1,:)',mask1(2,:)',ones(size(mask1,2),1),max(mask1(1,:)),gene.stop-gene.start+1);
        end;
        mm_reads = size(reads1, 1);

        mask1 = get_reads(bam_unf, gene.chr, gene.start, gene.stop, '0', 0, 1000, 50000, 2, 10);
        if isempty(mask1),
            reads1 = zeros(0,gene.stop-gene.start+1);
        else
            reads1 = sparse(mask1(1,:)',mask1(2,:)',ones(size(mask1,2),1),max(mask1(1,:)),gene.stop-gene.start+1);
        end;
        unf_reads = size(reads1, 1);

        genes(i).multi_frac = zeros(1, length(genes(i).transcripts));
        genes(i).multi_frac(:) = deal(mm_reads / unf_reads);
    end;

    save(filename, 'genes');

