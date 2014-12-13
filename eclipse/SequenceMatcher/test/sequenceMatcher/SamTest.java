package sequenceMatcher;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;

public class SamTest {

	@Test
	public void canConvertfastaListToSQHeader() {
		
		ArrayList<String[]> fastaList= new ArrayList<String[]>();
		fastaList.add(new String[] {"seq1", "ACTG"});
		fastaList.add(new String[] {"seq2", "ACTGACTG"});
		fastaList.add(new String[] {"seq3", "ACTGACTGACTG"});

		String expect= "@SQ	SN:seq1	LN:4\n"
   					 + "@SQ	SN:seq2	LN:8\n"
					 + "@SQ	SN:seq3	LN:12";

		String sqHeader= Sam.fastaListToSQHeader(fastaList);
		assertEquals(expect, sqHeader);
	}

	@Test
	public void canGetAlnStartPosFromRefAndReadStrings(){
		int s;
		String aln= "ACTG";
		String ref= "ACTG";
		s= Sam.getAlnStartPos(aln, ref);
		assertTrue(1 == s);
		
		aln= "NACTG";
		ref= "-ACTG";
		s= Sam.getAlnStartPos(aln, ref);
		assertTrue(1 == s);
		
		aln= "--TG";
		ref= "ACTG";
		s= Sam.getAlnStartPos(aln, ref);
		assertTrue(3 == s);
		
		aln= "---G";
		ref= "ACTG";
		s= Sam.getAlnStartPos(aln, ref);
		assertTrue(4 == s);

		aln= "ACTG";
		ref= "----";
		s= Sam.getAlnStartPos(aln, ref);
		assertTrue(-1 == s);

		aln= "----";
		ref= "----";
		s= Sam.getAlnStartPos(aln, ref);
		assertTrue(-1 == s);
		
	}
	
	@Test
	public void canGetCigarFromAln(){
		String aln;
		String ref;
		String cigar;

		aln= "ACTG";
		ref= "ACTG";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("4M", cigar);
		
		aln= "NACTG";
		ref= "-ACTG";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("1I4M", cigar);

		aln= "NNNACTG";
		ref= "---ACTG";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("3I4M", cigar);

		aln= "ACT-G";
		ref= "ACTNG";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("3M1D1M", cigar);

		aln= "ACT---G";
		ref= "ACTNNNG";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("3M3D1M", cigar);

		aln= "ACNTG";
		ref= "AC-TG";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("2M1I2M", cigar);

		aln= "ACNNNTG";
		ref= "AC---TG";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("2M3I2M", cigar);

		aln= "ACTGNN";
		ref= "ACTG--";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("4M2I", cigar);

		aln= "-C-G-N";
		ref= "A-T-N-";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("1I1D1I1D1I", cigar);

		aln= "--ACTGN--";
		ref= "NNNNACTGN";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("5M", cigar);

		aln= "N-CTTAGTTCAGGAGAT-GT--GAC----N";
		ref= "-NATAA---AGCTGATCGTTGAACCGTGGN";
		cigar= Sam.getCigarFromAln(aln, ref);
		assertEquals("1I1D4M3I8M1D2M2D3M4D1M", cigar);	
	}
	
	@Test
	public void canGetOriginalReadSeqFromAln(){
		String aln= "-ACTG--ACT";
		String read= Sam.getSeqFromAln(aln);
		assertEquals("ACTGACT", read);
	}
	
	@Test
	public void canConvertMatchToTags(){
		Match m= new Match();
		m.setSeqA("ACTGN");
		m.setSeqB("ACTGA");
		m.setStrand("+");
		m.setNameA("read1");
		m.setNameB("ref1");
		m.align();
		m.setHD();
		m.setLD();
		m.setJWD();

		ArrayList<String> tags= Sam.matchToTagList(m);
		
		assertTrue(tags.contains("NM:i:1"));
		assertTrue(tags.contains("AS:i:18"));
		assertTrue(tags.contains("XL:i:1"));
		assertTrue(tags.contains("XH:i:1"));		
		assertTrue(tags.contains("XJ:f:0.9200"));
		assertTrue(tags.contains("XP:f:0.8000"));
		assertTrue(tags.contains("XR:Z:ACTGN"));
		assertTrue(tags.contains("XS:Z:ACTGA"));
		
		Sam sam= Sam.matchToSam(m);
		
		System.out.println(sam.toString());
		
	}
}
