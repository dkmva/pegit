import regex

from design.helpers import reverse_complement, gc, dgn_to_regex


class Anzalone:

    @staticmethod
    def find_spacers(nuclease, reference_sequence, altered_sequence, start, end, spacer_search_range, cloning_method, **options):
        """Find candidate spacers for pegRNA selection.

        Finds all spacers with a cut site within spacer_search_range of the edit.
        Sorts spacers according to pam disruption, distance to edit and score.
        """
        spacers = []
        scoring_spacers = []
        sense = reference_sequence[:start + nuclease.downstream_from_cut_site]
        sense_offset = max(0, start - spacer_search_range - nuclease.cut_site_position)
        nucleotide_difference = len(altered_sequence) - len(reference_sequence)

        antisense = reverse_complement(
            reference_sequence[end - nuclease.downstream_from_cut_site - max(0, nucleotide_difference):])
        antisense_offset = end - max(0, nucleotide_difference) + min(
            len(antisense) - (spacer_search_range + nuclease.cut_site_position + nuclease.downstream_from_cut_site), 0)

        pam_motif = dgn_to_regex(nuclease.pam_motif) + '$'
        for match in regex.finditer(nuclease.target_motif,
                                    sense[-spacer_search_range - nuclease.cut_site_position - nuclease.downstream_from_cut_site:],
                                    regex.IGNORECASE, overlapped=True):
            spacer = match.group('spacer')
            pam = match.group('PAM')
            pos = sense_offset + match.start() + len(match.group('upstream'))
            pam_disrupted = not regex.search(pam_motif,
                                             altered_sequence[pos + len(spacer):pos + len(spacer) + len(pam)],
                                             regex.IGNORECASE)
            cut_site = pos + nuclease.cut_site_position - len(match.group('upstream'))
            distance = start - cut_site
            if nuclease.get_cloning_strategy(cloning_method).can_express(spacer):
                spacers.append({'spacer': spacer,
                                'position': pos,
                                'cut_site': cut_site,
                                'strand': 1,
                                'pam': (pam, pos + len(spacer)),
                                'pam_disrupted': pam_disrupted,
                                'distance': distance,
                                })

                scoring_spacers.append(match.group().upper())

        for match in regex.finditer(nuclease.target_motif, antisense[
                                                      -spacer_search_range - nuclease.cut_site_position - nuclease.downstream_from_cut_site:],
                                    regex.IGNORECASE, overlapped=True):
            spacer = match.group('spacer')
            pam = match.group('PAM')
            pos = antisense_offset + spacer_search_range - match.start() - len(
                match.group('upstream')) - 1 + nuclease.cut_site_position
            pam_disrupted = not regex.search(pam_motif, reverse_complement(altered_sequence[pos - len(spacer) - len(
                pam) + 1 + nucleotide_difference:pos + 1 + nucleotide_difference - len(spacer)]), regex.IGNORECASE)
            cut_site = pos - nuclease.cut_site_position + len(match.group('upstream')) + 1
            distance = cut_site - end + max(0, nucleotide_difference)
            if nuclease.get_cloning_strategy(cloning_method).can_express(spacer):
                spacers.append({'spacer': spacer,
                                'position': pos,
                                'cut_site': cut_site,
                                'strand': -1,
                                'pam': (pam, pos - len(spacer) - len(pam) + 1 + nucleotide_difference),
                                'pam_disrupted': pam_disrupted,
                                'distance': distance,
                                })

                scoring_spacers.append(match.group().upper())

        for i, score in enumerate(nuclease.score_spacers(scoring_spacers)):
            spacers[i]['score'] = score

        return sorted(spacers, key=lambda x: (not x['pam_disrupted'], x['distance'], x['score']))

    @staticmethod
    def find_nicking_spacers(nuclease, reference_sequence, altered_sequence, spacer_strand, cut_site, scaffold,
                             nicking_range, cloning_method, **options):
        """Find spacers for nicking the opposite strand."""
        spacers = []
        scoring_spacers = []

        nt_difference = len(altered_sequence) - len(reference_sequence)

        if spacer_strand == 1:
            reference_sequence = reverse_complement(reference_sequence)
            altered_sequence = reverse_complement(altered_sequence)
            cut_site = len(altered_sequence) - cut_site
        sequence = altered_sequence[
                   max(0, cut_site - nuclease.cut_site_position - nicking_range):
                   min(len(altered_sequence), cut_site + nuclease.downstream_from_cut_site + nicking_range)
                   ].upper()
        ref = reference_sequence[
              max(0, cut_site - nuclease.cut_site_position - nicking_range):
              min(len(reference_sequence), cut_site + nuclease.downstream_from_cut_site + nicking_range - nt_difference)
              ].upper()

        cut_site = nuclease.cut_site_position + nicking_range

        for match in regex.finditer(nuclease.target_motif, sequence, regex.IGNORECASE, overlapped=True):
            spacer = match.group('spacer')
            pos = match.start() + len(match.group('upstream'))
            wt_pos = pos
            cut = match.start() + nuclease.cut_site_position
            nick_location = cut_site - cut
            if nick_location < 0:
                wt_pos -= nt_difference

            kind = '3'
            wt_score = 1
            alt_bind = sequence[pos:pos + len(spacer) + len(nuclease.pam_motif)].upper()
            wt_bind = ref[wt_pos:wt_pos + len(spacer) + len(nuclease.pam_motif)].upper()
            if nuclease._is3b(alt_bind, wt_bind):
                kind = '3b'
                wt_score = nuclease._calc_wt_score(alt_bind, wt_bind)

            info = cloning_method.make_nicking_oligos(spacer, scaffold)
            info['position'] = nick_location
            info['spacer'] = spacer
            info['kind'] = kind
            info['wt_score'] = wt_score
            info['offset'] = nuclease.cut_site_position - len(match.group('upstream'))
            if nuclease.get_cloning_strategy(cloning_method).can_express(spacer):
                spacers.append(info)
                scoring_spacers.append(match.group().upper())

        for i, score in enumerate(nuclease.score_spacers(scoring_spacers)):
            spacers[i]['score'] = score

        return sorted(spacers, key=lambda x: (x['wt_score'], not (abs(x['position']) > 50), -x['score']))

    @staticmethod
    def _make_pbs_sequence(nuclease, reference, pbs_min_length, pbs_max_length, strategy, **options):
        """Find a suggested PBS length, and generate all possible PBS candidate lengths.

        Selects the shortest PBS sequence with a GC content in the range [0.4,0.6].
        If no sequence is within this range, selects the shortest PBS with a GC content closest to 0.5.
        """
        pbs_length = pbs_min_length - 1
        lengths = []
        while pbs_length < pbs_max_length:
            pbs_length += 1
            pbs = reference[-pbs_length:]
            if not nuclease.get_cloning_strategy(strategy).can_express(reverse_complement(pbs)):
                continue
            if 0.4 <= gc(pbs) <= 0.6:
                break
            lengths.append((abs(0.5 - gc(pbs)), len(pbs), pbs))
        else:
            try:
                pbs = sorted(lengths, key=lambda x: x[:1])[0][2]
            except IndexError:
                pbs = reference[-pbs_min_length:]

        # Create all possible PBS sequences within range limits.
        alt_lengths = [reference[-pbs_length:] for pbs_length in range(pbs_min_length, pbs_max_length + 1)]
        alt_lengths = [seq for seq in alt_lengths if nuclease.get_cloning_strategy(strategy).can_express(reverse_complement(seq))]
        return pbs, alt_lengths

    @staticmethod
    def _make_rt_sequence(nuclease, reference, cut_dist, nucleotide_difference, alteration_length, rt_min_length,
                          rt_max_length, strategy, **options):
        """Find a suggested RT template length, and generate alterniative RT template lengths."""
        rt_template_length = rt_min_length
        to_position = cut_dist + rt_template_length + nucleotide_difference
        rt_template = reference[:to_position]
        last_valid = rt_template
        # For large alterations, longer template is probably preferred
        while (nuclease._filter_extension(rt_template,
                                     strategy) or rt_template_length <= alteration_length * 2) and rt_template_length <= rt_max_length:
            rt_template += reference[to_position]
            if not nuclease._filter_extension(rt_template, strategy) and nuclease.get_cloning_strategy(strategy).can_express(reverse_complement(rt_template)):
                last_valid = rt_template
            to_position += 1
            rt_template_length += 1

        rt_template = last_valid
        # Create all possible RT templates within range limits.
        lengths = []
        for rt_template_length in range(rt_min_length, rt_max_length + 1):
            template = reference[:cut_dist + rt_template_length + nucleotide_difference]
            if not nuclease._filter_extension(template, strategy) and nuclease.get_cloning_strategy(strategy).can_express(reverse_complement(template)):
                lengths.append(template)
        return rt_template, lengths

    @classmethod
    def make_extension_sequence(cls, nuclease, reference_sequence, altered_sequence, spacer_strand, spacer_cut_site, cut_dist,
                                alteration_length, pbs_min_length, pbs_max_length, rt_min_length, rt_max_length,
                                strategy,
                                **options):
        """Create the pegRNA extension sequence.

        pegRNA extension sequences consist of a PBS and a RT template.
        The PBS is upstream of the cut site. The RT template is downstream of the cut site and contains the edit sequence.

        """
        nucleotide_difference = len(altered_sequence) - len(reference_sequence)
        if spacer_strand == 1:
            nucleotide_difference = min(0, nucleotide_difference)
            pbs_reference = reference_sequence[:spacer_cut_site]
            rt_reference = altered_sequence[spacer_cut_site:]
        else:
            pbs_reference = reverse_complement(reference_sequence[spacer_cut_site:])
            rt_reference = reverse_complement(altered_sequence[:spacer_cut_site + nucleotide_difference])
        pbs, pbs_lengths = cls._make_pbs_sequence(nuclease, pbs_reference.upper(), pbs_min_length, pbs_max_length, strategy,
                                                  **options)
        rt, rt_lengths = cls._make_rt_sequence(nuclease, rt_reference, cut_dist, nucleotide_difference, alteration_length,
                                               rt_min_length, rt_max_length, strategy, **options)
        pbs_length = len(pbs)
        rt_length = len(rt)

        # Generate all combinations of PBS and RT template sequences that are not identical to the primary suggestion.
        alternate_extensions = []
        for alt_pbs in pbs_lengths:
            alt_pbs_length = len(alt_pbs)
            alt_pbs_gc = round(gc(alt_pbs), 2)
            for alt_rt in rt_lengths:
                alt_rt_length = len(alt_rt)
                alt_rt_gc = round(gc(alt_rt), 2)
                if alt_pbs_length == pbs_length and alt_rt_length == rt_length:
                    continue
                alternate_extensions.append({'pbs_length': alt_pbs_length, 'rt_template_length': alt_rt_length,
                                             'sequence': reverse_complement(alt_pbs + alt_rt), 'pbs_gc': alt_pbs_gc,
                                             'rt_gc': alt_rt_gc})

        return pbs_length, rt_length, reverse_complement(pbs + rt), alternate_extensions