package bisReadBias;

import java.util.ArrayList;
import java.util.List;

import org.junit.*;

import com.google.common.collect.HashMultiset;

import static org.junit.Assert.*;

public class ReadProfileTest {
	
	@Test
	public void ReadProfileIsIncremented(){
		
		ReadProfile prof= new ReadProfile();
		prof.addReadMethyl("..M..u..**".getBytes(), true);
		prof.addReadMethyl("..M..M..**".getBytes(), true);
		
		assertEquals(2, prof.getRead1().get(2).count('M'));
		assertEquals(1, prof.getRead1().get(5).count('M'));
		
		// Read 2
		prof.addReadMethyl("..M..u..**".getBytes(), false);
		prof.addReadMethyl("..M..M..***".getBytes(), false);
		
		assertEquals(1, prof.getRead2().get(10).count('*'));
	}
	
	@Test
	public void acceptsAlignedReadObj() {
		AlignedRead aln= new AlignedRead();
		aln.setMethylbases("..M..u..**".getBytes());
		aln.setFirst(true);
		
		ReadProfile prof= new ReadProfile();
		prof.addReadMethyl(aln.getMethylbases(), aln.isFirst());
		assertEquals(1, prof.getRead1().get(2).count('M'));
	}
	
	@Test
	public void testPrintingToString(){
		
		ReadProfile prof= new ReadProfile();
		prof.addReadMethyl("..M..u..**".getBytes(), true);
		prof.addReadMethyl("..M..M..**".getBytes(), true);
		prof.addReadMethyl("..M..u..**".getBytes(), false);
		prof.addReadMethyl("..M..M..**".getBytes(), false);
		
		String printed= prof.toString();
		System.out.println(printed);		
	}
}
