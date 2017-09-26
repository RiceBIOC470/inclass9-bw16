% Inclass assignment 9

% The accession number for human NOTCH1 mRNA is AF308602
% 1. Read the information from this entry into matlab

accession='AF308602';
gb_data=getgenbank(accession);


%%
% 2. Write code that runs a blast query on the first 500 base pairs of this
% gene against the refseq_rna database

seq_begin=gb_data.Sequence(1:500);
[requestID,requestTime]=blastncbi(seq_begin,'blastn','Database','refseq_rna');
blast_data = getblast(requestID,'WaitTime',requestTime)

%%
% 3. Find the three highest scoring hits from other species and identify
% the length of the alignment and fraction of matches/mismatches. 

blast_data.Hits(1) %% Homo sapiens
blast_data.Hits(2) %% Pan troglodytes
blast_data.Hits(3) %% Pan troglodytes
blast_data.Hits(4) %% Pan troglodytes
blast_data.Hits(5) %% Pan troglodytes
blast_data.Hits(6) %% Rhinopithecus bieti
blast_data.Hits(7) %% Cercocebus atys

blast_data.Hits(2).HSPs.Alignment
blast_data.Hits(6).HSPs.Alignment
blast_data.Hits(7).HSPs.Alignment

Identity2=blast_data.Hits(2).HSPs(1).Identities
Identity6=blast_data.Hits(6).HSPs(1).Identities
Identity7=blast_data.Hits(7).HSPs(1).Identities

fraction2=Identity2/(500-Identity2)
fraction6=Identity6/(500-Identity6)
fraction7=Identity7/(500-Identity7)
%%
% 4. Run the same query against the database est_human. Comment on the
% sequences that you find. 
seq_begin=gb_data.Sequence(1:500);
[requestID2,requestTime2]=blastncbi(seq_begin,'blastn','Database','est_human');
blast_data_h=getblast(requestID2,'WaitTime',requestTime2)

blast_data_h.Hits(1) 
blast_data_h.Hits(2) 
blast_data_h.Hits(3)

blast_data_h.Hits(1).HSPs.Identities
blast_data_h.Hits(2).HSPs.Identities
blast_data_h.Hits(3).HSPs.Identities

%%Bingyan Wu:all alignments are against human genome