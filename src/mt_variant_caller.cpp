#include "mt_variant_caller.h"

// MtVariantCaller implementation
MtVariantCaller::MtVariantCaller(const Config& config) : _config(config) {
    // load fasta
    reference = _config.reference_file;

    _get_calling_interval();
    _print_calling_interval();

    // keep the order of '_samples_id' as the same as input bamfiles
    _get_sample_id_from_bam();
}

MtVariantCaller::~MtVariantCaller() {
    // 析构函数的实现，如果不需要特殊操作，则为空
}

bool MtVariantCaller::run() {

    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    bool is_success = _caller_process();

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();
    std::cout << "[INFO] " + ct + ". Done for Variant calling, "
              << difftime(now, real_start_time) << " (CPU time: "
              << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC << ") seconds elapsed."
              << std::endl;

    return is_success;
}

void MtVariantCaller::_get_sample_id_from_bam() {
    time_t real_start_time = time(0);

    // Loading sample ID in BAM/CRMA files from RG tag.
    if (_config.filename_has_samplename)
        std::cout << "[INFO] load samples id from filename directly, becuase you set "
                     "--filename-has-samplename\n";

    std::string samplename, filename;
    size_t si;
    for (size_t i(0); i < _config.bam_files.size(); ++i) {

        if ((i+1) % 1000 == 0)
            std::cout << "[INFO] loading "   << i+1 << "/" << _config.bam_files.size()
                      << " alignment files." << std::endl;

        if (_config.filename_has_samplename) {
            filename = ngslib::remove_filename_extension(ngslib::basename(_config.bam_files[i]));
            si = filename.find('.');
            samplename = si > 0 && si != std::string::npos ? filename.substr(0, si) : filename;
        } else {
            // Get sampleID from BAM header, a bit time-consuming.
            ngslib::BamHeader bh(_config.bam_files[i]);
            samplename = bh.get_sample_name();
        }

        if (!samplename.empty()) {
            _samples_id.push_back(samplename);
        } else {
            throw std::invalid_argument("[MtVariantCaller::_load_sample_id_from_bam] " +
                                        _config.bam_files[i] + " sample ID not found.\n");
        }
    }

    // check samples duplication
    std::vector<std::string> duplicate_samples = ngslib::find_duplicates(_samples_id);
    if (!duplicate_samples.empty()) {
        std::cout << "[WARNING] Find " << duplicate_samples.size() << " duplicated samples within "
                  << "the input bamfiles: " + ngslib::join(duplicate_samples, ",") + "\n"
                  << std::endl;
    }

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "[INFO] " + ct + ". Done for loading all samples' id from alignment files, "
              << difftime(now, real_start_time) << " seconds elapsed.\n"
              << std::endl;

    return;
}

void MtVariantCaller::_get_calling_interval() {

    std::vector<GenomeRegion> regions;
    if (!_config.calling_regions.empty()) {
        std::vector<std::string> rg_v;
        ngslib::split(_config.calling_regions, rg_v, ",");

        for (size_t i(0); i < rg_v.size(); ++i) {
            if (rg_v[i].length() == 0) continue; // ignore empty string
            regions.push_back(_make_genome_region(rg_v[i]));
        }
    } else {
        // Call the whole genome
        int n = reference.nseq();
        for (size_t i(0); i < n; ++i) {
            std::string ref_id = reference.iseq_name(i);
            regions.push_back(GenomeRegion(ref_id, 1, reference.seq_length(ref_id)));
        }
    }

    _calling_intervals.clear();
    // split region into small pieces by chunk_size
    for (size_t i(0); i < regions.size(); ++i) {
        uint32_t total_length = regions[i].end - regions[i].start + 1;
        for (uint32_t j(0); j < total_length; j += _config.chunk_size) {
            uint32_t start = regions[i].start + j;
            uint32_t end = std::min(regions[i].end, start + _config.chunk_size - 1);
            _calling_intervals.push_back(GenomeRegion(regions[i].ref_id, start, end));
        }
    }

    return;
}

GenomeRegion MtVariantCaller::_make_genome_region(std::string gregion) {

    // Genome Region, 1-based
    std::string ref_id;
    uint32_t reg_start, reg_end;

    std::vector<std::string> gr;
    ngslib::split(gregion, gr, ":");     // get reference id
    ref_id = gr[0];

    if (gr.size() == 2) {                // 'start-end' or start
        std::vector<uint32_t> gs;
        ngslib::split(gr[1], gs, "-");   // get position coordinate
        reg_start = gs[0];
        reg_end = (gs.size() == 2) ? gs[1] : reference.seq_length(ref_id);
    } else {
        reg_start = 1;
        reg_end = reference.seq_length(ref_id);  // the whole ``ref_id`` length
    }

    if (reg_start > reg_end) {
        throw std::invalid_argument("[ERROR] start postion is larger than end position in "
                                    "-r/--regions " + gregion);
    }

    return GenomeRegion(ref_id, reg_start, reg_end);  // 1-based
}

void MtVariantCaller::_print_calling_interval() {

    std::string ref_id;
    uint32_t reg_start, reg_end;
    std::cout << "---- Calling Intervals ----\n";
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        std::cout << i+1 << " - " << _calling_intervals[i].ref_id << ":"
                  << _calling_intervals[i].start << "-" << _calling_intervals[i].end << "\n";
    }
    std::cout << "\n";
    return;
}

bool MtVariantCaller::_caller_process() {

    // Get filepath and stem name
    std::string _bname = ngslib::basename(_config.output_file);
    size_t si = _bname.find(".vcf");
    std::string stem_bn = (si > 0 && si != std::string::npos) ? _bname.substr(0, si) : _bname;

    std::string outdir = ngslib::dirname(_config.output_file);
    std::string cache_outdir = outdir + "/cache_" + stem_bn;
    ngslib::safe_mkdir(cache_outdir);

    // Variant calling
    // 以区间为单位进行变异检测, 每个区间里先按照样本调用多线程，然后合并样本，多线程遍历位点并行处理
    bool is_success = true;
    std::vector<std::string> sub_vcf_files;
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        // 按样本进行 pileup，在函数里按样本进行多线程，记录每个样本在每个位点上的 pileup 信息
        PosMapVector samples_posmapinfo_v;
        samples_posmapinfo_v.reserve(_samples_id.size() + 1);  // reserve the memory before push_back
        bool is_empty = _fetch_base_in_region(_calling_intervals[i], samples_posmapinfo_v);
        if (is_empty) {
            std::cerr << "[WARNING] No reads in region: "   << _calling_intervals[i].ref_id << ":"
                      << _calling_intervals[i].start << "-" << _calling_intervals[i].end    << "\n";
            continue;
        }


        //////////////////////////////////////////////
        // Call variants in parallel
        std::string regstr = _calling_intervals[i].ref_id + "_" +
                             std::to_string(_calling_intervals[i].start) + "_" +
                             std::to_string(_calling_intervals[i].end);
        std::string sub_vcf_fn = cache_outdir + "/" + stem_bn + "." + regstr + ".vcf.gz";
        sub_vcf_files.push_back(sub_vcf_fn);
    }

    // Merge VCF
    std::string header = vcf_header_define(_config.reference_file, _samples_id);
    merge_file_by_line(sub_vcf_files, _config.output_file, header, true);

    const tbx_conf_t bf_tbx_conf = {1, 1, 2, 0, '#', 0};  // {preset, seq col, beg col, end col, header-char, skip-line}
    if ((ngslib::suffix_name(_config.output_file) == ".gz") &&          // create index
        tbx_index_build(_config.output_file.c_str(), 0, &bf_tbx_conf))  // file suffix is ".tbi"
    {
        throw std::runtime_error("tbx_index_build failed: Is the file bgzip-compressed? "
                                 "Check this file: " + _config.output_file + "\n");
    }

    return is_success;
}

bool MtVariantCaller::_fetch_base_in_region(const GenomeRegion gr, PosMapVector &samples_posmapinfo_v) {
    // The expend size of region, 100bp is enough
    static const uint32_t REG_EXPEND_SIZE = 100;
    uint32_t exp_reg_start = gr.start > REG_EXPEND_SIZE ? gr.start - REG_EXPEND_SIZE : 1;
    uint32_t exp_reg_end   = gr.end + REG_EXPEND_SIZE;
    std::string exp_regstr = gr.ref_id + ":" + std::to_string(exp_reg_start) + "-" + std::to_string(exp_reg_end);

    std::string fa_seq = this->reference[gr.ref_id];     // use the whole sequence of ``ref_id`` for simply
    ThreadPool thread_pool(this->_config.thread_count);  // set multiple-thread

    std::vector<std::future<PosMap>> pos_map_results;
    // Loop all alignment files
    for(size_t i(0); i < this->_config.bam_files.size(); ++i) {
        pos_map_results.push_back(thread_pool.submit(call_variant_in_sample, this->_config.bam_files[i],                                      
                                                     std::cref(fa_seq), gr, std::cref(this->_config)));
    }

    samples_posmapinfo_v.clear(); // clear the vector before push_back
    bool is_empty = true;
    for (auto && p: pos_map_results) { // Run and make sure all processes are finished
        // 一般来说，只有当 valid() 返回 true 的时候才调用 get() 去获取结果，这也是 C++ 文档推荐的操作。
        if (p.valid()) {
            // get() 调用会改变其共享状态，不再可用，也就是说 get() 只能被调用一次，多次调用会触发异常。
            PosMap pm = p.get();
            samples_posmapinfo_v.push_back(pm);

            if (!pm.empty()) is_empty = false;
        }
    }

    if (samples_posmapinfo_v.size() != this->_config.bam_files.size())
        throw std::runtime_error("[_fetch_base_in_region] 'samples_posmapinfo_v.size()' "
                                 "should be the same as '_config.bam_files.size()'");
    return is_empty;  // no cover reads in 'GenomeRegion' if empty.
}

// Seek the base information of each sample in the region.
PosMap call_variant_in_sample(const std::string sample_bam_fn, 
                              const std::string &fa_seq,
                              const GenomeRegion gr,
                              const MtVariantCaller::Config &config) {
    // The expend size of region, 100bp is enough.
    static const uint32_t REG_EXTEND_SIZE = 100;
    uint32_t extend_start     = gr.start > REG_EXTEND_SIZE ? gr.start - REG_EXTEND_SIZE : 1;
    uint32_t extend_end       = gr.end + REG_EXTEND_SIZE;
    std::string extend_regstr = gr.ref_id + ":" + std::to_string(extend_start) + "-" + std::to_string(extend_end);

    // 位点信息存入该变量, 且由于是按区间读取比对数据，key 值无需再包含 ref_id，因为已经不言自明
    PosMap sample_posinfo_map;           // key: position, value: alignment information
    ngslib::Bam bf(sample_bam_fn, "r");  // open bamfile in reading mode (one sample, one bamfile)
    if (bf.fetch(extend_regstr)) {       // Set 'bf' only fetch alignment reads in 'exp_regstr'.
        hts_pos_t map_ref_start, map_ref_end;  // hts_pos_t is uint64_t
        std::vector<ngslib::BamRecord> sample_target_reads; 
        ngslib::BamRecord al;  // alignment read
        while (bf.next(al) >= 0) {  // -1 => hit the end of alignement file.
            if (al.mapq() < config.min_mapq || al.is_duplicate() || al.is_qc_fail() ||
                (al.is_paired() && config.proper_pairs_only && !al.is_proper_pair()))
            {
                continue;
            }

            if (al.is_paired() && config.pairs_map_only) {
                std::string tid_name(al.tid_name(bf.header()));
                std::string mate_tid_name(al.mate_tid_name(bf.header()));

                // only use the paired reads which mapped to the same chromosome
                if (tid_name != mate_tid_name) continue;
            }

            map_ref_start = al.map_ref_start_pos() + 1;  // al.map_ref_start_pos() is 0-based, convert to 1-based
            map_ref_end   = al.map_ref_end_pos();        // al.map_ref_end_pos() is 1-based

            // Only fetch reads which in [reg_start, reg_end]
            if (gr.start > map_ref_end) continue;
            if (gr.end < map_ref_start) break;

            sample_target_reads.push_back(al);  // record all the proper reads of sample
        }

        if (sample_target_reads.size() > 0) {
            // get alignment information of [i] sample.
            seek_position(fa_seq, sample_target_reads, gr, sample_posinfo_map);
        }
    }

    // 信息抽提：先计算并返回每个样本在该区间里的突变信息，若不如此处理，内存吃不消
    // 这样做的好处还可为多样本 joint-calling 奠下可能. 
    for (auto &pos_align_info: sample_posinfo_map) {
        VariantInfo vi = variant_caller_unit(pos_align_info.second, config.heteroplasmy_threshold);
    }

    return sample_posinfo_map;
}

void seek_position(const std::string &fa_seq,   // must be the whole chromosome sequence
                   const std::vector<ngslib::BamRecord> &sample_map_reads,  // record the alignment reads of sample
                   const GenomeRegion gr,
                   PosMap &sample_posinfo_map)  // key: position, value: alignment information
{
    if (!sample_posinfo_map.empty())
        throw std::runtime_error("[seek_position] 'sample_posinfo_map' must be empty.");

    // A vector of: (cigar_op, read position, reference position, read base, read_qual, reference base)
    std::vector<ngslib::ReadAlignedPair> aligned_pairs;
    for (auto &al: sample_map_reads) {  // loop mapping reads
        AlignBase ab;
        ab.map_strand = al.map_strand();  // '*', '-' or '+'
        ab.mapq = al.mapq();

        // std::string ref_bases;
        uint32_t map_ref_pos;
        aligned_pairs = al.get_aligned_pairs(fa_seq);
        for (size_t i(0); i < aligned_pairs.size(); ++i) {
            map_ref_pos = aligned_pairs[i].ref_pos + 1;  // ref_pos is 0-based, convert to 1-based;

            if (gr.start > map_ref_pos) continue;
            if (gr.end < map_ref_pos) break;

            // 'BAM_XXX' are macros for CIGAR, which defined in 'htslib/sam.h'
            if (aligned_pairs[i].op == BAM_CMATCH ||  /* CIGAR: M */ 
                aligned_pairs[i].op == BAM_CEQUAL ||  /* CIGAR: = */
                aligned_pairs[i].op == BAM_CDIFF)     /* CIGAR: X */
            {
                // SNV. Only one character
                ab.ref_base  = aligned_pairs[i].ref_base[0];
                ab.read_base = aligned_pairs[i].read_base[0];
                ab.base_qual = aligned_pairs[i].read_qual[0];
            } else if (aligned_pairs[i].op == BAM_CINS) {  /* CIGAR: I */
                // Insertion
                if (!aligned_pairs[i].ref_base.empty()) {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got reference base in insertion region.");
                }

                // roll back one position to the leftmost to insertion break point.
                --map_ref_pos;
                ab.ref_base  = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base
                ab.read_base = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].read_base;

                // mean quality of the whole insertion sequence
                double total_score = 0;
                for (size_t i = 0; i < aligned_pairs[i].read_base.size(); ++i)
                    total_score += aligned_pairs[i].read_qual[i];
                ab.base_qual = uint8_t(total_score / aligned_pairs[i].read_base.size());

            } else if (aligned_pairs[i].op == BAM_CDEL) {  /* CIGAR: D */
                // Deletion.
                if (!aligned_pairs[i].read_base.empty()) {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got read bases in deletion region.");
                }

                // roll back one position to the leftmost to deletion break point.
                --map_ref_pos;
                ab.ref_base  = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].ref_base;
                ab.read_base = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base
                // set to be mean quality of the whole read if deletion
                ab.base_qual = uint8_t(al.mean_qqual()) + 33; // 33 is the offset of base QUAL; 

            } else { 
                // Skip other kinds of CIGAR symbals.
                continue;
            }

            // qpos is 0-based, conver to 1-based to set the rank of base on read.
            ab.rpr = aligned_pairs[i].qpos + 1;

            // 以 map_ref_pos 为 key，将所有的 read_bases 信息存入 map 中，多个突变共享同个 ref_pos，
            if (sample_posinfo_map.find(map_ref_pos) == sample_posinfo_map.end()) {
                // First level. If the position is not in the map, insert it.
                AlignInfo pos_align_info(gr.ref_id, map_ref_pos);
                sample_posinfo_map.insert({map_ref_pos, pos_align_info});
            }
            
            // collected all the base informations of the read into the map.
            sample_posinfo_map[map_ref_pos].align_bases.push_back(ab);
        }
    }

    return;
}

VariantInfo variant_caller_unit(const AlignInfo &pos_align_info, double min_af) {

    BatchInfo smp_bi;
    smp_bi.ref_id  = pos_align_info.ref_id;
    smp_bi.ref_pos = pos_align_info.ref_pos;
    for (auto &ab: pos_align_info.align_bases) {
        smp_bi.ref_bases.push_back(ab.ref_base);
        smp_bi.mapqs.push_back(ab.mapq);
        smp_bi.map_strands.push_back(ab.map_strand);

        smp_bi.align_bases.push_back(ab.read_base);
        smp_bi.align_base_quals.push_back(ab.base_qual);
        smp_bi.base_pos_ranks.push_back(ab.rpr);
    }

    std::transform(smp_bi.align_bases.begin(), smp_bi.align_bases.end(), smp_bi.align_bases.begin(), ::toupper);
    BaseType bt(&smp_bi, min_af);
    bt.lrt();  // use likelihood ratio test to detect candidate variant

    return get_variant(bt, &smp_bi);
}

VariantInfo get_variant(const BaseType &bt, const BatchInfo *smp_bi) {
    VariantInfo vi;

    return vi;
}
