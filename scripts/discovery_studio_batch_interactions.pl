#!/usr/bin/perl -w

use strict;
use MdmDiscoveryScript;

my ($manifest_path, $output_path, $status_path) = @ARGV;
die "Usage: discovery_studio_batch_interactions.pl MANIFEST.tsv OUTPUT.csv STATUS.csv\n"
    unless defined $status_path;

sub csv_value {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/\r?\n/ /g;
    $value =~ s/"/""/g;
    return '"' . $value . '"';
}

sub csv_row {
    my ($handle, @values) = @_;
    print $handle join(',', map { csv_value($_) } @values) . "\n";
}

sub bool_value {
    my ($value) = @_;
    return $value ? 'True' : 'False';
}

sub site_atoms {
    my ($site) = @_;
    return join(';', map { $_->FullName } @{$site->Atoms});
}

sub site_role {
    my ($site, $ligand) = @_;
    my $atoms = $site->Atoms;
    return '' unless $atoms->Count;
    return $atoms->Item(0)->IsMemberOf($ligand) ? 'ligand' : 'receptor';
}

open(my $manifest, '<', $manifest_path) or die "Could not open $manifest_path: $!\n";
open(my $output, '>', $output_path) or die "Could not write $output_path: $!\n";
open(my $status, '>', $status_path) or die "Could not write $status_path: $!\n";

csv_row($output, qw(
    structure_id role candidate result_folder interaction_type interaction_name
    interaction_categories category_names distance_A from_atoms to_atoms
    from_role to_role from_chemistry to_chemistry from_is_donor to_is_donor
    from_is_acceptor to_is_acceptor from_is_positive to_is_positive
    from_is_negative to_is_negative from_is_pi_system to_is_pi_system
));
csv_row($status, qw(structure_id role candidate status error interaction_count));

my $header = <$manifest>;
while (my $line = <$manifest>) {
    chomp $line;
    my ($structure_id, $role, $candidate, $result_folder, $pdb_path) = split(/\t/, $line, 5);
    my $count = 0;
    my $error = '';
    eval {
        my $document = DiscoveryScript::Open($pdb_path);
        die "Could not open staged PDB" unless $document;
        $document->FocusOnFirstLigand();
        my $ligand = $document->LigandInFocus;
        my $receptor = $document->ReceptorInFocus;
        die "No ligand was recognized" unless $ligand;
        die "No receptor was recognized" unless $receptor;

        my $monitor = $document->CreateLigandNonbondMonitor({
            IsMonitoringIntraligandNonbonds => False,
            RequirePredefinedLigand => True,
            AddHydrogens => False,
        });
        die "Could not create ligand nonbond monitor" unless $monitor;
        $document->UpdateObjects();

        foreach my $nonbond (@{$monitor->Nonbonds}) {
            my $from = $nonbond->FromSite;
            my $to = $nonbond->ToSite;
            my @types = @{$nonbond->NonbondTypes};
            my @categories = @{$nonbond->NonbondCategories};
            foreach my $type (@types) {
                csv_row($output,
                    $structure_id, $role, $candidate, $result_folder,
                    $type, $nonbond->TypeNames, join(';', @categories),
                    $nonbond->CategoryNames, $nonbond->Distance,
                    site_atoms($from), site_atoms($to),
                    site_role($from, $ligand), site_role($to, $ligand),
                    $from->ChemistryName, $to->ChemistryName,
                    bool_value($from->IsHydrogenDonor), bool_value($to->IsHydrogenDonor),
                    bool_value($from->IsHydrogenAcceptor), bool_value($to->IsHydrogenAcceptor),
                    bool_value($from->IsPositive), bool_value($to->IsPositive),
                    bool_value($from->IsNegative), bool_value($to->IsNegative),
                    bool_value($from->IsPiSystem), bool_value($to->IsPiSystem),
                );
                ++$count;
            }
        }
        eval { $document->Close(); };
    };
    if ($@) {
        $error = $@;
        $error =~ s/\s+$//;
        csv_row($status, $structure_id, $role, $candidate, 'failed', $error, $count);
    } else {
        csv_row($status, $structure_id, $role, $candidate, 'complete', '', $count);
    }
}

close($manifest);
close($output);
close($status);
