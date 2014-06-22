package bisReadBias;

import static org.junit.Assert.*;

import java.io.File;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.junit.*;
import org.mockito.Mockito;

import static org.mockito.Mockito.*;


public class AlignedReadTest {

	
	private SAMRecord rec= mock(SAMRecord.class);
	private AlignedRead alignedRead= new AlignedRead();
	
	@Test
	public void readBaseseAreCorrectlyReturned() throws Exception {

//		IndexedFastaSequenceFile faidx=new IndexedFastaSequenceFile(new File("test_data/chr7.fa"));
		SAMFileReader samfile= new SAMFileReader(new File("test_data/read01.sam"));
		SAMRecord rec= samfile.iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		
		AlignedRead alignedRead= new AlignedRead();
		String basesFromCigar= null;

		// Full match:
		rec.setCigarString("10M");
		rec.setReadBases("AAAAAAAAAA".getBytes());	
		basesFromCigar= new String(alignedRead.getAlignedBasesForRead(rec));
		assertEquals("AAAAAAAAAA", basesFromCigar);
		
		// With soft clip
		rec.setCigarString("2S5M3S");
		rec.setReadBases("AAAAAAAAAA".getBytes());
		basesFromCigar= new String(alignedRead.getAlignedBasesForRead(rec));
		assertEquals("NNAAAAANNN", basesFromCigar);

		//With skipped bases
		rec.setCigarString("7M10N3M");
		rec.setReadBases("AAAAAAACCC".getBytes());
		basesFromCigar= new String(alignedRead.getAlignedBasesForRead(rec));
		assertEquals("AAAAAAACCC", basesFromCigar);
		
		//With insertion to the reference
		rec.setCigarString("7M3I");
		rec.setReadBases("AAAAAAACCC".getBytes());
		basesFromCigar= new String(alignedRead.getAlignedBasesForRead(rec));
		assertEquals("AAAAAAACCC", basesFromCigar);

		//With deletion to the reference
		rec.setCigarString("10M3D");
		rec.setReadBases("AAAAAAACCC".getBytes());
		basesFromCigar= new String(alignedRead.getAlignedBasesForRead(rec));
		assertEquals("AAAAAAACCC", basesFromCigar);

		//Only clipped bases
		rec.setCigarString("10S");
		rec.setReadBases("AAAAAAACCC".getBytes());
		basesFromCigar= new String(alignedRead.getAlignedBasesForRead(rec));
		assertEquals("NNNNNNNNNN", basesFromCigar);
	
	}
	
	@Test
	public void referenceBaseseAreCorrectlyReturned() throws Exception {
		
		//AACCGGTTAANNAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
		IndexedFastaSequenceFile fa=new IndexedFastaSequenceFile(new File("test_data/chr7.small.fa"));		
		AlignedRead alignedRead= new AlignedRead();
		String basesFromRef= null;
		
		// ---------------------------------------------------------------------
		SAMRecord rec= new SAMFileReader(new File("test_data/read01.sam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		basesFromRef= new String(alignedRead.getReferenceBasesForRead(rec, fa));
		assertEquals("AACCGGTTAA", basesFromRef);
		
		// ---------------------------------------------------------------------
		rec= new SAMFileReader(new File("test_data/read02.sam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		basesFromRef= new String(alignedRead.getReferenceBasesForRead(rec, fa));
		assertEquals("AACCGGTNNN", basesFromRef);
		
		// ---------------------------------------------------------------------
		rec= new SAMFileReader(new File("test_data/read03.sam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		basesFromRef= new String(alignedRead.getReferenceBasesForRead(rec, fa));
		assertEquals("NNNGGTTAAN", basesFromRef);
		
		// ---------------------------------------------------------------------
		rec= new SAMFileReader(new File("test_data/read04.sam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		basesFromRef= new String(alignedRead.getReferenceBasesForRead(rec, fa));
		assertEquals("GGNNAACCGG", basesFromRef);
	}
	
	@Test
	public void readMappedToRevStrandEqualsFastqSeq() throws Exception{
		
		AlignedRead alignedRead= new AlignedRead();
		String bases= null;
		
		// ---------------------------------------------------------------------
		SAMRecord rec= new SAMFileReader(new File("test_data/read01.rev.bam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		
		AlignedRead aln= new AlignedRead(rec);
		bases= new String(aln.getReadbases());
		assertEquals("CAATGGACAA", bases);		
	
	}
	
	@Test
	public void mate2MappedToForwStrandEqualsFastqSeq() throws Exception{
		
		AlignedRead alignedRead= new AlignedRead();
		String bases= null;
		
		// ---------------------------------------------------------------------
		SAMRecord rec= new SAMFileReader(new File("test_data/read01.R2.forw.sam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		
		AlignedRead aln= new AlignedRead(rec);
		bases= new String(aln.getReadbases());
		assertEquals("CATTGTCCAT", bases);
	}
	
	@Test
	public void mate2MappedToRevStrandEqualsFastqSeq() throws Exception{
		
		AlignedRead alignedRead= new AlignedRead();
		String bases= null;
		
		// ---------------------------------------------------------------------
		SAMRecord rec= new SAMFileReader(new File("test_data/read01.R2.rev.sam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		
		AlignedRead aln= new AlignedRead(rec);
		bases= new String(aln.getReadbases());
		
		assertEquals("GCAAATCGAG", bases);
	}
	
	@Test
	public void referenceForRead1RevStrandEqualsFastqSeq() throws Exception{
		
		IndexedFastaSequenceFile fa=new IndexedFastaSequenceFile(new File("test_data/ref.fa"));
		AlignedRead alignedRead= new AlignedRead();
		String bases= null;
		
		// ---------------------------------------------------------------------
		SAMRecord rec= new SAMFileReader(new File("test_data/read01.rev.bam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		
		AlignedRead aln= new AlignedRead(rec, fa);;
		assertArrayEquals("CAATGGACAA".getBytes(), aln.getRefbases());
	}
	
	@Test
	public void shouldReturnMethylString() throws Exception{
		
		byte[] read= "CAATGGAcAA".getBytes();
		byte[] ref=  "cAATGGACAA".getBytes();
		AlignedRead aln= new AlignedRead(); 
		byte[] bisStr= aln.methylArray(read, ref);
		assertArrayEquals("M......M..".getBytes(), bisStr);

		AlignedRead aln2= new AlignedRead();
		aln2.setReadbases("ACTG".getBytes());
		aln2.setRefbases( "ACTG".getBytes());
		assertArrayEquals(     ".M..".getBytes(), aln2.methylArray());
		
	}
	
	@Test
	public void toStringShouldReturnALine() throws Exception{
		IndexedFastaSequenceFile fa=new IndexedFastaSequenceFile(new File("test_data/ref.fa"));
		// ---------------------------------------------------------------------
		SAMRecord rec= new SAMFileReader(new File("test_data/read01.rev.bam")).iterator().next();
		rec.setValidationStringency(ValidationStringency.SILENT);
		rec.setFlags(89);
		AlignedRead aln= new AlignedRead(rec, fa);
		
		assertEquals("CAATGGACAA	CAATGGACAA	M......M..	false	true", aln.toString());
	}
}