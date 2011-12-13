SET default_tablespace = pg_default;
drop table pileup.cage20101217;
select read_table($$ file:'F:/data/20101217_CAGE050810/cage_050810_bwa_20101216.clean.single.forward.slice.pileup', 
    header:True, dest_table:'forward'  $$);
create table cage20101217 AS(select *, '+'::text AS strand from "forward");
drop table "forward";
alter table cage20101217 set tablespace hdd_free_agent;
alter table cage20101217 set schema pileup;

select read_table($$ file:'F:/data/20101217_CAGE050810/cage_050810_bwa_20101216.clean.single.reverse.slice.pileup', 
    header:True, dest_table:'reverse' $$);
insert into pileup.cage20101217 (select *, '-'  from reverse);
drop table "reverse";

create index ind_cage20101217 on pileup.cage20101217(pos, rname, strand);
alter index pileup.ind_cage20101217_forward set tablespace hdd_free_agent;
vacuum analyze pileup.cage20101217;

alter table pileup.cage20101217 alter column strand TYPE text;

select * from pileup.cage20101217_reverse where rname = '7' and pos between 27658519 and 27661277;

select * from pileup.cage20101217 inner join sus_scrofa_sscrofa9_59_gtf on
              pileup.cage20101217.rname = sus_scrofa_sscrofa9_59_gtf.rname AND
              pileup.cage20101217.strand = sus_scrofa_sscrofa9_59_gtf.strand AND
              pileup.cage20101217.pos between f_start and f_end
where gtf_attribute(attributes, 'gene_name') in ('TNFA', 'CSF1R')
order by sus_scrofa_sscrofa9_59_gtf.rname, f_start;

where gtf_attribute(attributes, 'gene_name') in ('TNFA_PIG')
, 'CSF1R'

, '%CSF1R%'
select * from sus_scrofa_sscrofa9_59_gtf where gtf_attribute(attributes, 'gene_name') in ('TNFA_PIG');
" gene_id "ENSSSCG00000001404"; transcript_id "ENSSSCT00000001533"; exon_number "1"; gene_name "TNF"; transcript_name "TNFA_PIG";"
