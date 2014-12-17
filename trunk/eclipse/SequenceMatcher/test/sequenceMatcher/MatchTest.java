package sequenceMatcher;

import static org.junit.Assert.*;

import org.junit.Test;

public class MatchTest {

	@Test
	public void canGetDistance() {
		//Create a set of sequences and a single sequence. 
		String[] seqA= {"seq1-A", "ACTGA"};
		String[] seqB= {"seq1-B", "ACTGN"};

		Match matches= new Match(seqA, seqB);	
		
		matches.computeHD(0);
		assertTrue(matches.getHD() == -1);
		
		matches.computeLD(1);
		assertTrue(matches.getLD() == 1);
		
		matches.computeJWD();
		assertTrue(matches.getJWD() == 0.92);
		
		double d= matches.getFilterDistance("HD", -1);
		assertTrue(d == 1);

		d= matches.getFilterDistance("LD", 1);
		assertTrue(d == 1);
		
		d= matches.getFilterDistance("LD", -1);
		assertTrue(d == 1);
		
		d= matches.getFilterDistance("LD", 0);
		assertTrue(d == -1);
		
		matches= new Match(seqA, seqA);
		d= matches.getFilterDistance("LD", 0);
		assertTrue(d == 0);

		d= matches.getFilterDistance("LD", 1);
		assertTrue(d == 0);

		d= matches.getFilterDistance("LD", -1);
		assertTrue(d == 0);
		
	}
	
	@Test
	public void canConvertStringToMatchObject(){
		//seq_A	seq_B strand LD HD JWD len_A len_B NM aln_score len_aln n_ident n_sim pct_ident aln_A aln_B
		String line= "seq1	seq2	-	16	20	0.74	30	30	22	2	14	16	0.39	36	NT------CGTGCGGCTATACGTGCAGGGCATTCAT	ATGAATGCCCTGCACGTATAGCCGCACG------AN";
		Match m= new Match(); 
		m.stringToMatchObj(line);
		
//		System.out.println(m);
		
		assertEquals("seq1", m.getNameA());
		assertEquals("seq2", m.getNameB());
		assertEquals("-", m.getStrand());
		assertEquals(16, m.getLD());
		assertEquals(20, m.getHD());
		assertEquals(0.74, m.getJWD(), 0.001);
		assertEquals(30, m.getLen_A());
		assertEquals(30, m.getLen_B());
		assertEquals(22, m.getNM());
//		// ...
		assertEquals("NT------CGTGCGGCTATACGTGCAGGGCATTCAT", m.getAlnA());
		assertEquals("ATGAATGCCCTGCACGTATAGCCGCACG------AN", m.getAlnB());
	}
	
	@Test
	public void canConvertMatchToString(){
		Match m= new Match();
		// System.out.println(m.toString());
		m.setNameA("myseq");
		// System.out.println(m.toString());
		assertTrue(m.toString().startsWith("myseq"));
	}
	
	@Test
	public void testLineToMatchObjAndBack(){
	
		String line= "seq1	seq2	-	16	20	0.74	30	30	22	2	14	16	0.39	36	NT------CGTGCGGCTATACGTGCAGGGCATTCAT	ATGAATGCCCTGCACGTATAGCCGCACG------AN";

		Match m= new Match(); 
		m.stringToMatchObj(line);
		
		String newLine= m.toString();
		assertEquals(line, newLine);
	}
}
