create temp table mart_tome_h37 (
  ensembl_gene_id varchar,
  ensembl_transcript_id varchar,
  chr varchar,
  strand varchar,
  transcript_start int,
  transcript_end int
  )

copy mart_tome_h37 from 'C:/Downloads/mart_export_h37.txt' with csv header delimiter E',';

create temp table mart_tome_ss9 (
  ensembl_gene_id varchar,
  ensembl_transcript_id varchar,
  chr varchar,
  strand varchar,
  transcript_start int,
  transcript_end int
  )

copy mart_tome_ss9 from 'C:/Downloads/mart_export_ss9.txt' with csv header delimiter E',';

select * from mart_tome_ss9 limit 10;

drop table if exists mart_tome;
create table mart_tome AS (
  select 'h37' AS genome, mart_tome_h37.* from mart_tome_h37
  union select 'ss9' AS genome, mart_tome_ss9.* from mart_tome_ss9
  );

alter table mart_tome add constraint PK_mart_tome_genome_transcript primary key (genome, ensembl_transcript_id);

-- Reset strand from varchar to boolean
update mart_tome set strand = true
  where strand like '1'; 
update mart_tome set strand = false
  where strand like '-1'; 

update mart_tome set strand = strand::boolean;

alter table mart_tome rename column strand to strand_old;

alter table mart_tome add column strand boolean;

update mart_tome set strand = strand_old::boolean;

alter table mart_tome drop column strand_old;
-------------------------------------------------------------------------------
comment on table mart_tome is 'Ensembl transcripts imported from Biomart';
comment on column mart_tome.genome is 'Genome from which transcriptome has been retrieved. ss9 = S.scrofa 9, h37 = H. sapiens (GRCh37)'

drop table if exists mart_tome_h37;
drop table if exists mart_tome_ss9;


