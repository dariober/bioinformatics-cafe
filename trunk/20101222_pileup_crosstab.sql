set work_mem= '512MB';

-- drop table pileup_rnaseq_lps;
create table pileup_rnaseq_lps AS(
  select * from pileup_rnaseq_ctrl where 1=2
  );
alter table pileup_rnaseq_lps set tablespace hdd_free_agent;
comment on table pileup_rnaseq_lps is 'RNAseq from 30/06/2010 in pileup format, LPS library only';
copy pileup_rnaseq_lps from 'D:/Tritume/20100630_RNAseq_LPS.pileup' with csv header delimiter E'\t';
create index ind_rnaseq_lps on pileup_rnaseq_lps(pos, rname);
alter index ind_rnaseq_lps set tablespace hdd_free_agent;
-- vacuum analyze pileup_rnaseq_lps;



---------------------------[ Generate cross-tab of pileups ]-------------------

-- drop table pileup_rnaseq;
create table pileup_rnaseq AS(
  select rname, pos, ref_base, NULL::int AS n_ctrl, 
                               NULL::int AS n_lps from pileup_rnaseq_ctrl
  union 
  select rname, pos, ref_base, NULL::int AS n_ctrl, 
                               NULL::int AS n_lps from pileup_rnaseq_lps
  );
create index ind_rnaseq_all on pileup_rnaseq(rname, pos) tablespace hdd_free_agent;
vacuum analyze pileup_rnaseq;

update pileup_rnaseq set n_ctrl = read_count from pileup_rnaseq_ctrl 
                                             where pileup_rnaseq.rname = pileup_rnaseq_ctrl.rname and	
                                                   pileup_rnaseq.pos = pileup_rnaseq_ctrl.pos; 
update pileup_rnaseq set n_lps = read_count from pileup_rnaseq_lps 
                                             where pileup_rnaseq.rname = pileup_rnaseq_lps.rname and
                                                   pileup_rnaseq.pos = pileup_rnaseq_lps.pos;
