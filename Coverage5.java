// import sys
// import re
// import pysam
// from sys import argv # argv[1]: prefix, argv[2]: PATH of BAM file
// from os import path
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.ValidationStringency;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.BlockingQueue;

public class Coverage5 {
    SAMSequenceDictionary seqDict = null;
    static String rg = null;
    public void process_chunk(SAMRecord linechunk){
    }

    public static void main(String[] argv){
        String prefix = "";
        String bamfile = argv[0];
        int concurrency = Integer.parseInt(argv[1]);
        if(argv.length > 2){
            rg = argv[2];
        }
        Coverage5 cov = new Coverage5();
        try {
            final SamReader reader = SamReaderFactory.makeDefault().open(new File(bamfile));
            SAMFileHeader header = reader.getFileHeader();
            cov.seqDict = header.getSequenceDictionary();
            for(int i = 0; i<cov.seqDict.size(); i++){
                SAMSequenceRecord seqRec = cov.seqDict.getSequence(i);
                System.out.println(seqRec.getSequenceName() + " " + seqRec.getSequenceLength() + " " + seqRec.getSequenceIndex());
            }
            if(reader.hasIndex()){
                System.err.println("Index is OK");
            }else {
                System.err.println("No index");
                System.exit(0);
            }
            reader.close();
            cov.run(bamfile, concurrency);
        }catch(Exception e){
            e.printStackTrace();
        }
    }   
    public void run(String bamfile, int concurrency){
        LinkedBlockingQueue<Coverage> queue = new LinkedBlockingQueue<Coverage>(concurrency);
        for(int i = 0; i<seqDict.size(); i++){
/*
            if(!seqDict.getSequence(i).getSequenceName().equals("chrY")){
                continue;
            }*/
            Coverage c = new Coverage(bamfile, queue, seqDict.getSequence(i));
            try {
                queue.put(c);
                c.start();
            }catch(InterruptedException e){
                e.printStackTrace();
            }
        }
    }
    public class Coverage extends Thread{
        BlockingQueue queue;
        File bam;
        int[] depth_dic;
        SAMSequenceRecord srec;
        public Coverage(String bamfile, BlockingQueue q, SAMSequenceRecord srec_){
            queue = q;
            bam = new File(bamfile);
            srec = srec_;
        }
        /*
        linechunk = []
        prev_pos = -1
        prev_chrom = ""
        
        depth_dic = {}
        
        def process_chunk(linechunk):
            for chrom, pos, cigar in linechunk:
                cur_pos = pos
                pos_back = 0
                leftflag = True
                for s, n in cigar:
                    if not leftflag or (s != 4 and s != 5):
                        for pos_diff in range(n):
                            if (s != 1 and s != 4 and s != 5):
                                if (cur_pos-pos_back) in depth_dic:
                                    depth_dic[cur_pos-pos_back] += 1
                                else:
                                    depth_dic[cur_pos-pos_back] = 1
                            if s == 1:
                                pos_back += 1
                            cur_pos += 1
                        leftflag = False
        */
        public void process_chunk(ArrayList<Chunk> linechunk){
            for(Chunk chunk: linechunk){
                int cur_pos = chunk.pos;
                int pos_back = 0;
                boolean leftflag = true;

                for(CigarElement c: chunk.cigar){
                    if(!leftflag || !c.getOperator().isClipping()){
                        for(int pos_diff = 0; pos_diff < c.getLength(); pos_diff++){
                            if (c.getOperator() != CigarOperator.I && !c.getOperator().isClipping()){
                                depth_dic[cur_pos - pos_back]++;
                            }
                            if(c.getOperator() == CigarOperator.I) {
                                pos_back++;
                            }
                            cur_pos++;
                        }
                        leftflag = false;
                    }
                }
            }
        }
        public void run(){
            System.err.println("running " + srec.getSequenceName());
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bam);
            SAMFileHeader header = reader.getFileHeader();
            if(rg != null){
            	header.addReadGroup(new htsjdk.samtools.SAMReadGroupRecord("1"));
            }
            SAMRecordIterator it = reader.query(srec.getSequenceName(), 1, srec.getSequenceLength(), true);
            try {
                PrintWriter out = new PrintWriter(new java.io.BufferedWriter(new java.io.FileWriter(srec.getSequenceName() + "_depth.txt")));
                int prev_pos = -1;
                depth_dic = new int[srec.getSequenceLength() + 1];
                ArrayList<Chunk> linechunk = new ArrayList<Chunk>();
                while(it.hasNext()){
                    SAMRecord sam = it.next();
                    int pos = sam.getAlignmentStart();
                    Cigar cigar = sam.getCigar();
                    if(pos != prev_pos){
                        if(prev_pos != -1){
                            process_chunk(linechunk);
                            for(int i = prev_pos; i<pos; i++){
                                out.printf("%s:%d\t%d\n", srec.getSequenceName(), i, depth_dic[i]);
                                depth_dic[i] = 0;
                            }
                        }
                        prev_pos = pos;
                        linechunk.clear();
                        linechunk.add(new Chunk(srec.getSequenceName(), pos, cigar));
                    }else{
                        linechunk.add(new Chunk(srec.getSequenceName(), pos, cigar));
                    }
                }
                out.close();
                reader.close();
            }catch(IOException e){
                e.printStackTrace();
            }
            finish();
        }
        public void finish(){
            queue.remove(this);
        }
        public class Chunk {
            String chrom;
            int pos;
            Cigar cigar;
            public Chunk(String chr, int p, Cigar cgr){
                chrom = chr;
                pos = p;
                cigar = cgr;
            }
        }
    }
/*               
    prefix = ""
    if len(argv) > 2:
        prefix = argv[1]
        fn = argv[2]
    else:
        fn = argv[1]

    f = pysam.Samfile(fn, "rb")
    depth_dic = {}
    for cnt, ar in enumerate(f):
        if ar.tid == -1:
            continue
        chrom = f.getrname(ar.tid)
        if prev_chrom != chrom:
            if prev_chrom != "":
                fo.close()
        
            prev_chrom = chrom
            if prefix == "":
                fo = open(chrom+"_depth.txt", "w")
            else:
                fo = open(prefix+chrom+"_depth.txt", "w")
        pos = ar.pos+1
        cigar = ar.cigar
        
        if pos != prev_pos:
            if prev_pos != -1:
                process_chunk(linechunk)
                for i in range(prev_pos, pos):
                    try:
                        fo.write('%s:%d\t%d\n'%(chrom, i, depth_dic[i]))
                        del depth_dic[i]
                    except KeyError:
                        fo.write('%s:%d\t%d\n'%(chrom, i, 0))
            prev_pos = pos
            linechunk = [ (chrom, pos, cigar) ]
        else:
            linechunk.append( (chrom, pos, cigar) )
        if cnt % 50000 == 0:
            msg = "@"+chrom+": "+str(pos)
            sys.stdout.write(msg);sys.stdout.flush()
            sys.stdout.write("\b"*len(msg));sys.stdout.flush()
        

    process_chunk(linechunk)
    for i in range(pos, pos+150):
        try:
            fo.write('%s:%d\t%d\n'%(chrom, i, depth_dic[i]))
        except KeyError:
            fo.write('%s:%d\t%d\n'%(chrom, i, 0))

    f.close()
    fo.close()    
*/
}
