import uproot
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Define colors and markers for different PF candidate types, GenParticles, LostTracks, and Jets
pf_particle_colors = {
    211: ('red', 'o'),       # Charged hadron
    11: ('blue', 's'),       # Electron
    13: ('green', '^'),      # Muon
    22: ('yellow', 'v'),     # Photon
    130: ('orange', 'D'),    # Neutral hadron
    1: ('purple', 'x'),      # HF hadron
    2: ('cyan', '*'),        # HF EM particle
}

gen_particle_colors = {
    211: ('magenta', 'o'),   # Charged hadron
    11: ('lime', 's'),       # Electron
    13: ('darkgreen', '^'),  # Muon
    22: ('gold', 'v'),       # Photon
    130: ('darkorange', 'D'),# K0_L
    310: ('blue', 'D'),      # K0_S
    311: ('cyan', 'D'),      # K0
    111: ('brown', 'd'),     # Pi0
    1000015: ('black', 'P'), # SUSY stau (PDG ID 1000015)
    15: ('pink', '8'),       # Tau
}

lost_track_color = ('gray', 'h') # Lost tracks color and marker
jet_color = ('blue', '--')  # Jet color and linestyle

def determine_tau_decay_mode(types_gen, tau_children_mask):

    n_charged_hadrons = sum(1 for i, is_child in enumerate(tau_children_mask) if is_child and abs(types_gen[i]) == 211)
    n_photons = sum(1 for i, is_child in enumerate(tau_children_mask) if is_child and abs(types_gen[i]) == 22)
    n_electrons = sum(1 for i, is_child in enumerate(tau_children_mask) if is_child and abs(types_gen[i]) == 11)
    n_muons = sum(1 for i, is_child in enumerate(tau_children_mask) if is_child and abs(types_gen[i]) == 13)
    n_neutral_hadrons = sum(1 for i, is_child in enumerate(tau_children_mask) if is_child and abs(types_gen[i]) in [111, 130, 310, 311])

    if n_muons > 0:
        return "Muon Decay"
    elif n_electrons > 0:
        return "Electron Decay"
    elif n_charged_hadrons == 1 and n_photons == 0 and n_neutral_hadrons == 0:
        return "1-Prong Hadronic"
    elif n_charged_hadrons == 1 and n_neutral_hadrons > 0:
        return "1-Prong + Neutral Hadrons"
    elif n_charged_hadrons == 3 and n_neutral_hadrons == 0:
        return "3-Prong Hadronic"
    elif n_charged_hadrons == 3 and n_neutral_hadrons > 0:
        return "3-Prong + Neutral Hadrons"
    else:
        return "Other"


def plot_eta_phi(eta_pf, phi_pf, types_pf, energies_pf,
                 eta_gen, phi_gen, types_gen, energies_gen,
                 eta_lost, phi_lost, tau_index, event,
                 jet_eta=None, jet_phi=None, tau_children_mask=None, tau_vertex_x=None,
                 tau_vertex_y=None, tau_vertex_z=None, tau_eta=None, tau_phi=None):

    print(f"Plotting event {event}, tau {tau_index}")
    print(tau_eta,jet_eta)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(36, 18))

    # Plot PF candidates and LostTracks
    unique_labels_pf = set()
    # Dictionary mapping PDG IDs to their names
    pdg_names = {
        211: 'Charged Hadron',
        11: 'Electron',
        13: 'Muon',
        22: 'Photon',
        130: 'Neutral Hadron',
        310: 'Neutral Hadron',
        311: 'Neutral Hadron',
        111: 'Neutral Hadron',
        1000015: 'SUSY Stau',
        15: 'Tau'
    }
    
    for i in range(len(eta_pf)):
        particle_type = abs(types_pf[i])
        color, marker = pf_particle_colors.get(particle_type, ('black', 'o'))  # Default to black and circle marker if undefined type
        pdg_name = pdg_names.get(particle_type, 'Unknown')
        label = f'{pdg_name} ({particle_type})' if particle_type not in unique_labels_pf else ""
        ax1.scatter(phi_pf[i], eta_pf[i], color=color, marker=marker, alpha=0.7, s=300, label=label)
        unique_labels_pf.add(particle_type)

    # Plot LostTracks if available
    if len(eta_lost) > 0 and len(phi_lost) > 0:
        ax1.scatter(phi_lost, eta_lost, color=lost_track_color[0], marker=lost_track_color[1], alpha=0.5, s=300, label='LostTrack')

    # Plot reconstructed jet as a circle contour if available
    if jet_eta is not None and jet_phi is not None:
        circle = plt.Circle((jet_phi, jet_eta), 0.4, color=jet_color[0], fill=False, linestyle=jet_color[1], linewidth=2, label='Jet (R=0.4)')
        ax1.add_patch(circle)

    # Function to normalize phi difference
    def normalize_phi(phi_diff):
        return (phi_diff + np.pi) % (2 * np.pi) - np.pi

    # Function to calculate delta phi
    def delta_phi(phi1, phi2):
        dphi = phi1 - phi2
        while dphi > np.pi:
            dphi -= 2 * np.pi
        while dphi < -np.pi:
            dphi += 2 * np.pi
        return dphi

    # Highlight the 4 most energetic PF candidates with arrows and their energy
    if len(energies_pf) > 0:
        most_energetic_pf_indices = np.argsort(energies_pf)[-4:]
        directions = [(1, 1), (-1, 1), (-1, -1), (1, -1)]  # Default directions to avoid overlap
        for idx, i in enumerate(most_energetic_pf_indices):
            if jet_eta is not None and jet_phi is not None:
                # Calculate direction vector towards the jet
                dx = - delta_phi(jet_phi, phi_pf[i])
                dy = - jet_eta + eta_pf[i]
                norm = np.sqrt(dx**2 + dy**2)
                # if norm 0, set dx and dy to 1
                if norm == 0:
                    dx = 0.2
                    dy = 0.2
                else:
                    dx /= norm
                    dy /= norm
                    dx *= 0.2  # Fixed length for the arrow
                    dy *= 0.2
            # else:
            #     # Use default direction
            #     dx, dy = directions[idx]
            #     dx *= 0.1
            #     dy *= 0.1
                ax1.annotate(f'{energies_pf[i]:.2f} GeV', xy=(phi_pf[i], eta_pf[i]), xytext=(phi_pf[i] + dx, eta_pf[i] + dy),
                            arrowprops=dict(facecolor='blue', shrink=0.05), fontsize=23, color='blue')

    ax1.set_title(f"Reconstructed, Event {event}", fontsize=36)
    ax1.set_xlabel("Phi", fontsize=34)
    ax1.set_ylabel("Eta", fontsize=34)
    ax1.tick_params(axis='both', which='major', labelsize=28)
    ax1.grid(True)

    # Add legend for PF candidates, LostTracks, and Jets
    ax1.legend(fontsize=30)

    # Plot GenParticles
    unique_labels_gen = set()
    plotted_stau = False
    for i in range(len(eta_gen)):
        particle_type = abs(types_gen[i])
        if particle_type in {12, 14, 16}:  # Skip neutrinos
            continue
        if particle_type == 1000015 and plotted_stau:  # Skip duplicate stau particles
            continue
        if particle_type == 1000015:
            plotted_stau = True
        color, marker = gen_particle_colors.get(particle_type, ('black', 'o'))  # Default to black and circle marker if undefined type
        pdg_name = pdg_names.get(particle_type, 'Unknown')
        label = f'{pdg_name} ({particle_type})' if particle_type not in unique_labels_gen else ""
        if tau_children_mask is not None and tau_children_mask[i]:
            ax2.scatter(phi_gen[i], eta_gen[i], color=color, marker=marker, alpha=0.7, s=300, label=label, edgecolors='black', linewidths=2)
        else:
            ax2.scatter(phi_gen[i], eta_gen[i], color=color, marker=marker, alpha=0.7, s=300, label=label)
        unique_labels_gen.add(particle_type)

    # Highlight the 4 most energetic tau children with arrows and their energy
    if len(energies_gen) > 0:
        most_energetic_gen_indices = [i for i in np.argsort(energies_gen)[-4:] if tau_children_mask[i]]
        for idx, i in enumerate(most_energetic_gen_indices):
            # Calculate direction vector towards the tau center
            dx = - delta_phi(tau_phi,phi_gen[i])
            dy = - eta_gen[i] + tau_eta
            norm = np.sqrt(dx**2 + dy**2)
            if norm == 0:
                dx = 0.2
                dy = 0.2
            else:
                dx /= norm
                dy /= norm
                dx *= 0.2  # Fixed length for the arrow
                dy *= 0.2
            ax2.annotate(f'{energies_gen[i]:.2f} GeV', xy=(phi_gen[i], eta_gen[i]), xytext=(phi_gen[i] + dx, eta_gen[i] + dy),
                         arrowprops=dict(facecolor='blue', shrink=0.05), fontsize=23, color='blue')

    ax2.set_title(f"Generator, Event {event}", fontsize=36)
    
    if tau_vertex_x is not None and tau_vertex_y is not None and tau_vertex_z is not None:
        ax2.text(0.95, 0.15, f"Tau Vertex (r_xy, z) [cm]: ({np.sqrt(tau_vertex_x**2 + tau_vertex_y**2):.2f}, {tau_vertex_z:.2f})", transform=ax2.transAxes, fontsize=28, verticalalignment='bottom', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))

    ax2.set_xlabel("Phi", fontsize=34)
    ax2.set_ylabel("Eta", fontsize=34)
    ax2.tick_params(axis='both', which='major', labelsize=28)
    ax2.grid(True)

    # Add legend for GenParticles
    ax2.legend(fontsize=30)

    # Determine tau decay mode and add text to plot
    if tau_children_mask is not None:
        tau_decay_mode = determine_tau_decay_mode(types_gen, tau_children_mask)
        ax2.text(0.95, 0.05, f"Decay Mode: {tau_decay_mode}", transform=ax2.transAxes, fontsize=28,
                 verticalalignment='bottom', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))

    # find gen particle with pdg 15 and status 2
    tau_energy = None
    for i in range(len(eta_gen)):
        if abs(types_gen[i]) == 15:
            tau_energy = energies_gen[i]
            break

    # Add text to plot with tau energy
    if tau_energy is not None:
        ax2.text(0.95, 0.10, f"Tau Energy: {tau_energy:.2f} GeV", transform=ax2.transAxes, fontsize=28,
                 verticalalignment='bottom', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))

    # Align eta and phi axes on both plots, with a bit larger limits to see particles on the borders
    eta_min = min(np.min(eta_pf), np.min(eta_gen), np.min(eta_lost) if len(eta_lost) > 0 else np.inf, jet_eta - 0.5 if jet_eta is not None else np.inf) - 0.2
    eta_max = max(np.max(eta_pf), np.max(eta_gen), np.max(eta_lost) if len(eta_lost) > 0 else -np.inf, jet_eta + 0.5 if jet_eta is not None else -np.inf) + 0.2
    phi_min = min(np.min(phi_pf), np.min(phi_gen), np.min(phi_lost) if len(phi_lost) > 0 else np.inf, jet_phi - 0.5 if jet_phi is not None else np.inf) - 0.2
    phi_max = max(np.max(phi_pf), np.max(phi_gen), np.max(phi_lost) if len(phi_lost) > 0 else -np.inf, jet_phi + 0.5 if jet_phi is not None else -np.inf) + 0.2

    ax1.set_xlim(phi_min, phi_max)
    ax1.set_ylim(eta_min, eta_max)
    ax2.set_xlim(phi_min, phi_max)
    ax2.set_ylim(eta_min, eta_max)

    plt.tight_layout()
    if os.path.isdir("display") is False:
        os.makedirs("display")
    plt.savefig(f"display/event_{event}_tau_{tau_index}.png")
    plt.close()

def main():
    # Get the file path and number of events to analyze from command line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <path_to_root_file> <number_of_events>")
        sys.exit(1)

    file_path = sys.argv[1]
    num_events = int(sys.argv[2])

    # Open the ROOT file
    file = uproot.open(file_path)
    tree = file["Events"]

    # Define the branches to extract
    branches = [
        "GenPart_vertexX",
        "GenPart_vertexY",
        "GenPart_vertexZ",
        "PFCandidate_eta",
        "PFCandidate_phi",
        "PFCandidate_pdgId",
        "PFCandidate_pt",
        "PFCandidate_mass",  # Add this line
        "GenPart_eta",
        "GenPart_phi",
        "GenPart_pdgId",
        "GenPart_status",
        "GenPart_statusFlags",
        "GenPart_pt",
        "GenPart_mass",  # Add this line
        "GenVisTau_eta",
        "GenVisTau_phi",
        "LostTrack_eta",
        "LostTrack_phi",
        "Jet_eta",
        "Jet_phi",
    ]

    # Read the branches
    data = tree.arrays(branches, library="np")

    # Extract the variables
    eta_pf = data["PFCandidate_eta"]
    phi_pf = data["PFCandidate_phi"]
    pdg_id_pf = data["PFCandidate_pdgId"]
    pt_pf = data["PFCandidate_pt"]
    mass_pf = data["PFCandidate_mass"]  # Add this line

    eta_gen = data["GenPart_eta"]
    phi_gen = data["GenPart_phi"]
    pdg_id_gen = data["GenPart_pdgId"]
    status_gen = data["GenPart_status"]
    status_flags_gen = data["GenPart_statusFlags"]
    pt_gen = data["GenPart_pt"]
    mass_gen = data["GenPart_mass"]  # Add this line
    
    eta_vis_tau = data["GenVisTau_eta"]
    phi_vis_tau = data["GenVisTau_phi"]
    
    eta_lost = data["LostTrack_eta"]
    phi_lost = data["LostTrack_phi"]
    eta_jet = data["Jet_eta"]
    phi_jet = data["Jet_phi"]

    # Iterate over events and GenVisTaus and plot the PF candidates and GenParticles close to each
    for event_index in range(min(num_events, len(eta_vis_tau))):
        for tau_index in range(len(eta_vis_tau[event_index])):
            tau_eta = eta_vis_tau[event_index][tau_index]
            tau_phi = phi_vis_tau[event_index][tau_index]
            # Find GenPart corresponding to tau or stau (only last copy)
            for gen_index in range(len(eta_gen[event_index])):
                if abs(pdg_id_gen[event_index][gen_index]) == 15 and status_gen[event_index][gen_index] == 2 and (status_flags_gen[event_index][gen_index] & (1 << 12)):  # Check if GenVisTau is a tau or stau

                    # Define delta R condition for matching
                    delta_r_pf = np.sqrt((eta_pf[event_index] - tau_eta)**2 + (phi_pf[event_index] - tau_phi)**2)
                    delta_r_gen = np.sqrt((eta_gen[event_index] - tau_eta)**2 + (phi_gen[event_index] - tau_phi)**2)
                    delta_r_lost = np.sqrt((eta_lost[event_index] - tau_eta)**2 + (phi_lost[event_index] - tau_phi)**2)

                    # Mask for PF candidates, GenParticles, and LostTracks close to GenVisTau
                    mask_pf = delta_r_pf < 0.5
                    mask_gen = delta_r_gen < 0.5
                    mask_lost = delta_r_lost < 0.5

                    eta_pf_tau = eta_pf[event_index][mask_pf]
                    phi_pf_tau = phi_pf[event_index][mask_pf]
                    types_pf_tau = pdg_id_pf[event_index][mask_pf]
                    pt_pf_tau = pt_pf[event_index][mask_pf]
                    eta_gen_tau = eta_gen[event_index][mask_gen]
                    phi_gen_tau = phi_gen[event_index][mask_gen]
                    types_gen_tau = pdg_id_gen[event_index][mask_gen]
                    pt_gen_tau = pt_gen[event_index][mask_gen]
                    eta_lost_tau = eta_lost[event_index][mask_lost]
                    phi_lost_tau = phi_lost[event_index][mask_lost]
                    mass_pf_tau = mass_pf[event_index][mask_pf]
                    mass_gen_tau = mass_gen[event_index][mask_gen]

                    # Determine which GenParticles are children of tau
                    tau_children_mask = [(status_flags_gen[event_index][mask_gen][i] & (1 << 4)) != 0 for i in range(len(eta_gen_tau))]

                    # Find reconstructed jet that matches GenVisTau
                    jet_match_idx = None
                    for jet_index in range(len(eta_jet[event_index])):
                        delta_r_jet = np.sqrt((eta_jet[event_index][jet_index] - tau_eta)**2 + (phi_jet[event_index][jet_index] - tau_phi)**2)
                        if delta_r_jet < 0.3:
                            jet_match_idx = jet_index
                            break

                    jet_eta = eta_jet[event_index][jet_match_idx] if jet_match_idx is not None else None
                    jet_phi = phi_jet[event_index][jet_match_idx] if jet_match_idx is not None else None

                    # calculate energyies of pf candidates and gen particles
                    energies_pf = np.sqrt((pt_pf_tau*np.cosh(eta_pf_tau))**2 + mass_pf_tau**2)
                    energies_gen = np.sqrt((pt_gen_tau*np.cosh(eta_gen_tau))**2 + mass_gen_tau**2)

                    if jet_eta is not None and jet_phi is not None:
                        # Plot eta-phi distribution for PF candidates, LostTracks, GenParticles, and Jet close to the current GenVisTau or stau
                        plot_eta_phi(eta_pf_tau, phi_pf_tau, types_pf_tau, energies_pf, eta_gen_tau, phi_gen_tau, types_gen_tau, energies_gen, eta_lost_tau, phi_lost_tau, tau_index, event_index, jet_eta, jet_phi, tau_children_mask, data["GenPart_vertexX"][event_index][gen_index], data["GenPart_vertexY"][event_index][gen_index], data["GenPart_vertexZ"][event_index][gen_index], tau_eta, tau_phi)

if __name__ == "__main__":
    main()
