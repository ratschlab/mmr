function genes = add_transcript_length(fname),

    load(fname);

    for i = 1:length(genes),
        if mod(i, 100) == 0,
            fprintf(1, '.');
            if mod(i, 1000) == 0,
                fprintf(1, '%i\n', i);
            end;
        end;
        genes(i).transcript_length = [];
        for j = 1:length(genes(i).transcripts),
            if ~isempty(genes(i).exons{j}),
                genes(i).transcript_length(end + 1) = sum(genes(i).exons{j}(:, 2) - genes(i).exons{j}(:, 1) + 1);
            else
                genes(i).transcript_length(end + 1) = 0;
            end;
        end;
    end;

    save(fname, 'genes');
