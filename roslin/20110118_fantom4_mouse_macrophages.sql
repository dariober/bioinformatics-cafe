show search_path;
drop table "FANTOM4".mouse_macrophage_mapping;
select read_table($$ file: 'F:/data/FANTOM4/Tables/mouse/CAGE/mapping/con_mapping.tbl.txt', table:'FANTOM4.mouse_macrophage_mapping', header: True, limit: -1 $$);
select * from "FANTOM4".mouse_macrophage_mapping limit 1000;
COMMENT ON TABLE "FANTOM4".mouse_macrophage_mapping IS 'Mouse cage tags mapped from FANTOM4 see 18/01/2011'
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.id IS $$extracted tag sequence$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.library_count IS $$read count in this library. SUM( rna_count ) = library_count$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.edit_string IS $$required edition of the sequence for its matching to the genome$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.chrom IS $$the name of the chromosome$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.strand IS $$the strand in the chromosome$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.start IS $$the starting position of the feature in the chromosome$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.end IS $$the ending position of the feature in the chromosome$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.percentage IS $$percent match (always 100)$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.map_pos IS $$number of locations where this tag sequence mapped to the genome$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.ribo_flag IS $$tag mapped to ribosome equal or better than chromosome?$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.refseq_flag IS $$tag mapped to refseq better than chromosome?$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.rna IS $$RNA name$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.rna_count IS $$read count in the RNA$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.tpm_in_ribo IS $$tag per million per timecourse including ribosome flagged tags$$;
COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.tpm_ex_ribo IS $$tag per million per timecourse excluding ribosome_flagged tags$$;
-- COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.weight IS $$weight of expression. This is just internal value, and please ignore this column for  subsequent analysis$$;
-- COMMENT ON COLUMN "FANTOM4".mouse_macrophage_mapping.rescue_weight IS $$weight computed with the 'rescue'  strategy ( Faulkner et al. Genomics. 2008 91:281 )$$;

select distinct dataset_id from "FANTOM4".mouse_macrophage_mapping;

-- Select rows containing unstimulated samples:
create temp table mm_unstimulated AS(
  select * from "FANTOM4".mouse_macrophage_mapping 
  where rna in (
   'CBU macrophage unstimulated sample',
   'CBY macrophage cultured overnight in CSF-1(macrophage growth factor)',
   'CDQ macrophage cultured overnight in CSF-1(macrophage growth factor)',
   'CDV macrophage',
   'CGO macrophage unstimulated sample'
   )
);
select * from mm_unstimulated limit 10; -- 998347 rows

-- Produce TSSs:
drop table "FANTOM4".mouse_macrophage_tss_bed;
create table "FANTOM4".mouse_macrophage_tss_bed AS (
select 
    min( chrom ) AS chrom, 
    min( CASE WHEN strand like '+' THEN "start" 
         WHEN strand like '-' THEN "end" END ) AS tss_pos, -- Position of TSS 0-based coordinate.
    chrom || '_' || (CASE WHEN strand like '+' THEN "start" ELSE "end" END) || '_' || strand AS tss_id,
    min( strand ) AS strand, 
    rna,
    sum(rna_count)::int AS rna_count,
    sum(tpm_ex_ribo) AS tpm_ex_ribo,
    count(id)::int AS no_tags
FROM "FANTOM4".mouse_macrophage_mapping
GROUP BY rna, chrom || '_' || (CASE WHEN strand like '+' THEN "start" ELSE "end" END) || '_' || strand -- this is: GROUP BY tss_id
);
comment on table "FANTOM4".mouse_macrophage_tss_bed IS 'BED table of transcr. s. sites. Coordinates are 0-based. So, use tss_end if comparing to ensembl files. See 20110118_fantom4_mouse_macrophages.sql';

select * from "FANTOM4".mouse_macrophage_tss_bed limit 2000;

-------------------------------[ Tritume ]-------------------------------------

-- Produce TSSs:
drop table mm_tss;
create temp table mm_tss AS(
select rna, chrom, strand,
    CASE WHEN strand like '+' THEN "start" 
         WHEN strand like '-' THEN "end" END AS tss,
    sum(rna_count) as rna_count,
    sum(tpm_ex_ribo) as tpm,
    count(id) as no_ids
from mm_unstimulated
group by rna, chrom, strand, CASE WHEN strand like '+' THEN "start" WHEN strand like '-' THEN "end" END
); -- 639839 rows

select * from mm_tss order by no_ids desc limit 20;


select * from "FANTOM4".mouse_macrophage_mapping where 
  rna like 'CGO macrophage unstimulated sample' and chrom like 'chr11' and strand like '-' and ("end" = 108873394 or "start" = 108873394);


"CGO macrophage unstimulated sample";"chr11";"-";108873394;3696;0;75

