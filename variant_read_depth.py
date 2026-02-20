#!/usr/bin/env python3
import pysam
import os
import multiprocessing as mp


def process_variant(args):
    bam_path, chrom, pos, ref, alt = args

    bam = pysam.AlignmentFile(bam_path, "rb")

    hp1_ref = hp1_alt = 0
    hp2_ref = hp2_alt = 0
    unassigned = 0

    for col in bam.pileup(
        chrom,
        pos,
        pos + 1,
        truncate=True,
        stepper="samtools"
    ):
        if col.reference_pos != pos:
            continue

        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue

            aln = pr.alignment
            base = aln.query_sequence[pr.query_position]
            hp = aln.get_tag("HP") if aln.has_tag("HP") else None

            if base == ref:
                if hp == 1:
                    hp1_ref += 1
                elif hp == 2:
                    hp2_ref += 1
                else:
                    unassigned += 1
            elif base == alt:
                if hp == 1:
                    hp1_alt += 1
                elif hp == 2:
                    hp2_alt += 1
                else:
                    unassigned += 1

    bam.close()

    dp_hap = hp1_ref + hp1_alt + hp2_ref + hp2_alt + unassigned

    return hp1_ref, hp1_alt, hp2_ref, hp2_alt, unassigned, dp_hap


def define_depth(bam_path, vcf_in_path, vcf_out_path):

    vcf_in = pysam.VariantFile(vcf_in_path)

    # ---- header ----
    header = vcf_in.header.copy()
    header.info.add("HP1_REF", 1, "Integer", "Reads supporting REF on haplotype 1")
    header.info.add("HP1_ALT", 1, "Integer", "Reads supporting ALT on haplotype 1")
    header.info.add("HP2_REF", 1, "Integer", "Reads supporting REF on haplotype 2")
    header.info.add("HP2_ALT", 1, "Integer", "Reads supporting ALT on haplotype 2")
    header.info.add("UNASSIGNED", 1, "Integer", "Reads without HP tag")
    header.info.add("DP_HAP", 1, "Integer", "Total reads used for haplotype counting")

    vcf_out = pysam.VariantFile(vcf_out_path, "w", header=header)

    # ---- extract plain-python variant info ----
    variants = []
    records = []

    for rec in vcf_in:
        records.append(rec)

        if len(rec.alleles) != 2:
            variants.append(None)
        else:
            chrom = rec.contig
            pos = rec.pos - 1
            ref, alt = rec.alleles
            variants.append((bam_path, chrom, pos, ref, alt))

    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))

    with mp.Pool(nproc) as pool:
        results = pool.map(
            process_variant,
            [v for v in variants if v is not None]
        )

    res_iter = iter(results)

    # ---- write output ----
    for rec, var in zip(records, variants):

        if var is None:
            vcf_out.write(rec)
            continue

        hp1_ref, hp1_alt, hp2_ref, hp2_alt, unassigned, dp_hap = next(res_iter)

        new_rec = vcf_out.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=rec.alleles,
            id=rec.id,
            qual=rec.qual,
            filter=rec.filter,
            info=dict(rec.info)
        )

        new_rec.info["HP1_REF"] = hp1_ref
        new_rec.info["HP1_ALT"] = hp1_alt
        new_rec.info["HP2_REF"] = hp2_ref
        new_rec.info["HP2_ALT"] = hp2_alt
        new_rec.info["UNASSIGNED"] = unassigned
        new_rec.info["DP_HAP"] = dp_hap

        vcf_out.write(new_rec)

    vcf_in.close()
    vcf_out.close()


def main():
    vers = "2"
    samp = "wf_hvar"

    in_bam = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{vers}/{samp}_ROI/haplotagged.bam"
    in_vcf = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{vers}/{samp}_ROI/rn_{samp}.wf_snp.vcf.gz"
    out_vcf = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{vers}/{samp}_ROI/AD_{samp}.wf_snp.vcf.gz"

    define_depth(in_bam, in_vcf, out_vcf)


if __name__ == "__main__":
    main()

