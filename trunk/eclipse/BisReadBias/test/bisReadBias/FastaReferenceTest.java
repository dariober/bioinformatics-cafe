package bisReadBias;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.*;

public class FastaReferenceTest {

	@Test
	public void testPrinting(){
		FastaReference fa= new FastaReference(new File("test_data/ref.fa"));
		System.out.println(fa.toString());
	}
	
	@Test
	public void getSubsequenceFromFasta(){
		
		FastaReference fa= new FastaReference(new File("test_data/ref.fa"));
		byte[] subseq= fa.getSubsequenceAt("ref", 1, 10, false);		
		assertArrayEquals("TACGCACGGG".getBytes(), subseq);
		
		subseq= fa.getSubsequenceAt("ref", 2, 11, false);		
		assertArrayEquals("ACGCACGGGT".getBytes(), subseq);
	}

	@Test
	public void getFirstBases(){
		
		FastaReference fa= new FastaReference(new File("test_data/ref.fa"));
		byte[] subseq= fa.getSubsequenceAt("ref", 1, 3, false);
		assertArrayEquals("TAC".getBytes(), subseq);
		
	}

	@Test
	public void getLastBases(){
		
		FastaReference fa= new FastaReference(new File("test_data/ref.fa"));
		int refLength= fa.getReference().get("ref").length;
		byte[] subseq= fa.getSubsequenceAt("ref", 98, refLength, false);
		assertArrayEquals("TTC".getBytes(), subseq);
		
	}
	
	@Test
	public void getSubsequencePassedTheEndOfRef(){
		
		FastaReference fa= new FastaReference(new File("test_data/ref.fa"));
		byte[] subseq= fa.getSubsequenceAt("ref", 98, 104, true);
		
		assertArrayEquals("TTCNNNN".getBytes(), subseq);
	}
	
	
	
	@Test
	public void profileGetSubsequenceAtIsFastEnough(){
		
		FastaReference fa= new FastaReference(new File("test_data/ref.fa"));
		byte[] subseq= fa.getSubsequenceAt("ref", 1, 10, true);		
		assertArrayEquals("TACGCACGGG".getBytes(), subseq);
		
		long t0= System.currentTimeMillis();
		int i= 0;
		while(i < 10000000){
			subseq= fa.getSubsequenceAt("ref", 1, 10, true);
			i++;
		}
		long t1= System.currentTimeMillis();
//		System.out.println((t1-t0)/(float)1000);
		assertTrue((t1-t0) < 2000);
		
	}
}
