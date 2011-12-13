/*
  Custom annotation of affymetrix chip porcine na30.
  See LabBook 10 May 10
*/
-- drop table tmp_affy_transcriptome;
-- This file produced by 20100510_bowtie_porcine_na30.sh
select readbowtie('C:/Tritume/porcine_na30_Sus_scrofa.Sscrofa9.56.cdna.all.bowtiemap', 'tmp_affy_transcriptome');
select * from tmp_affy_transcriptome limit 20;

create temp table affy_tome AS(
-- Extract probe_set_id and probe_id from full fasta name
select substring(read_name 
                 from position('probe_set_id:' in read_name)+13 for 
                 position('|' in read_name) - (position('probe_set_id:' in read_name)+13)
                 ) AS probe_set_id,
       substring(read_name 
                 from position('probe_id:' in read_name)+9 for 
                 length(read_name) - (position('probe_id:' in read_name)+8)
                 ) AS probe_id,        
* from tmp_affy_transcriptome
);
select * from affy_tome limit 10;

-- Extract annotate probes not previously reported from Biomart
drop table new_affy;
create temp table new_affy AS(
select 
   affy_tome.probe_set_id, 
   affy_tome.refseq_name AS transcript_id, 
   count(probe_id) AS mapped_probes, -- Number of probes from probe_set_id matching this transcript
   no_probes -- Total number of probes in this probe_set_id
from affy_tome left join porcine_na30_ensembl on 
     -- Link to Biomart annotation
     affy_tome.probe_set_id = porcine_na30_ensembl.probe_set_id 
inner join 
     -- Number of probes in each probe_set_id
     (select count(probe_id) AS no_probes, probe_set_id from affy_tome group by probe_set_id) AS t1 on
     t1.probe_set_id = affy_tome.probe_set_id     
where porcine_na30_ensembl.probe_set_id is null
group by affy_tome.probe_set_id, affy_tome.refseq_name, no_probes
);

copy new_affy to 'C:/Tritume/affy_annotation.txt' with csv header delimiter E'\t';

-----------------------------------[ Tritume ]------------------------------------

select probe_set_id, count(transcript_id) from new_affy group by probe_set_id;
select distinct transcript_id from new_affy;


select count(distinct read_name) from tmp_affy_transcriptome limit 10;
select count(distinct refseq_name) from tmp_affy_transcriptome limit 10;

select read_name, count(distinct refseq_name) AS hits  from tmp_affy_transcriptome group by read_name having count(distinct refseq_name) > 1;

select * from tmp_affy_transcriptome where read_name like 'AFFX-Ss_Actin_5_s_at';