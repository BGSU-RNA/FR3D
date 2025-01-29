% Read a file of atom to atom mappings between standard nucleotides and modified nucleotides

function [modified_base_to_parent,modified_base_atom_list,modified_atom_to_parent,parent_atom_to_modified,modified_base_to_hydrogens,modified_base_to_hydrogen_coordinates] = zDefineModifiedNucleotides()

    % modified_base_to_parent = containers.Map();
    % modified_base_atom_list = containers.Map();
    % modified_atom_to_parent = containers.Map();
    % parent_atom_to_modified = containers.Map();
    % modified_base_to_hydrogens = containers.Map();
    % modified_base_to_hydrogen_coordinates = containers.Map();

    % load characteristics of standard bases
    zStandardBases;
    base_to_atoms = containers.Map();
    base_to_coordinates = containers.Map();
    base_to_atom_to_coordinates = containers.Map();
    % base atom names
    base_to_atoms('A') = A_Atoms;
    base_to_atoms('C') = C_Atoms;
    base_to_atoms('G') = G_Atoms;
    base_to_atoms('U') = U_Atoms;
    base_to_coordinates('A') = A_Stand;
    base_to_coordinates('C') = C_Stand;
    base_to_coordinates('G') = G_Stand;
    base_to_coordinates('U') = U_Stand;

    % store atom coordinates
    all_bases = base_to_atoms.keys();
    for i = 1:length(all_bases)
        base = all_bases{i};
        atom_to_coordinates = containers.Map();
        atoms = base_to_atoms(base);
        coordinates = base_to_coordinates(base);
        for k = 1:length(atoms)
            atom_to_coordinates(atoms{k}) = coordinates(k,:);
        end
        base_to_atom_to_coordinates(base) = atom_to_coordinates;
    end

    % Get current file path
    current_file = mfilename('fullpath');
    current_path = fileparts(current_file);
    filename = fullfile(current_path, 'atom_mappings.txt');

    % Read atom to atom mappings for modified nucleotides
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = lines{1};

    % store data by modified nucleotide
    modified_atom_map = containers.Map();
    for i = 1:length(lines)
        fields = strsplit(lines{i});
        if length(fields) == 4
            modified_nucleotide = fields{3};
            if ~isKey(modified_atom_map, modified_nucleotide)
                modified_atom_map(modified_nucleotide) = {};
            end
            modified_atom_map(modified_nucleotide) = [modified_atom_map(modified_nucleotide); {fields{1}, fields{2}, fields{4}}];
        end
    end

    % map DA, DC, DG, DT to A, C, G, U as if DNA nucleotides are modified RNA nucleotides
    modified_atom_map('DA') = modified_atom_map('0A');
    modified_atom_map('DC') = modified_atom_map('0C');
    modified_atom_map('DG') = modified_atom_map('0G');
    modified_atom_map('DT') = modified_atom_map('0U');

    % Initialize outputs
    modified_base_to_parent = containers.Map();
    modified_base_atom_list = containers.Map();
    modified_atom_to_parent = containers.Map();
    parent_atom_to_modified = containers.Map();
    modified_base_to_hydrogens = containers.Map();
    modified_base_to_hydrogen_coordinates = containers.Map();

    % Process the modified_atom_map
    modified_nucleotides = modified_atom_map.keys();
    for i = 1:length(modified_nucleotides)
        modified_nucleotide = modified_nucleotides{i};

        mam = modified_atom_map(modified_nucleotide);
        parent_nucleotide = mam{1,1};

        % fprintf('Processing %-3s with parent %2s\n', modified_nucleotide,parent_nucleotide);

        % deal with DNA bases by mapping them to RNA bases
        if ismember(parent_nucleotide, {'DA', 'DC', 'DG'})
            parent_nucleotide = parent_nucleotide(2);
        elseif parent_nucleotide == 'DT'
            parent_nucleotide = 'U';
        end

        modified_base_to_parent(modified_nucleotide) = mam{1,1};

        % make temporary lists for this modified base
        this_modified_base_atom_list = [];
        this_modified_atom_to_parent = containers.Map();
        this_parent_atom_to_modified = containers.Map();
        this_modified_base_to_hydrogens = [];
        this_modified_base_to_hydrogen_coordinates = containers.Map();
        parent_atom_to_hydrogen_coordinates = base_to_atom_to_coordinates(parent_nucleotide);

        % Process each atom mapping for the modified nucleotide
        for k = 1:size(mam, 1)
            parent_atom       = mam{k,2};
            modified_atom     = mam{k,3};

            if ismember(modified_nucleotide,{'DA','DC','DG','DT'})
                % replace HN with N in parent atom names
                if parent_atom == "O2'"
                    continue
                end
                modified_atom = strrep(modified_atom, 'HN', 'N');
            end

            % Map the modified atom to the parent atom
            this_modified_atom_to_parent(modified_atom) = parent_atom;
            this_parent_atom_to_modified(parent_atom) = modified_atom;

            % If parent atom is in the base
            if ismember(parent_atom, base_to_atoms(parent_nucleotide))
                this_modified_base_atom_list = [this_modified_base_atom_list, modified_atom];
                this_modified_atom_to_parent(modified_atom) = parent_atom;
                this_parent_atom_to_modified(parent_atom) = modified_atom;
                if parent_atom(1) == 'H'
                    this_modified_base_to_hydrogens = [this_modified_base_to_hydrogens, modified_atom];
                    this_modified_base_to_hydrogen_coordinates(modified_atom) = parent_atom_to_hydrogen_coordinates(parent_atom);
                end
            end
        end
        modified_atom_to_parent(modified_nucleotide) = this_modified_atom_to_parent;
        parent_atom_to_modified(modified_nucleotide) = this_parent_atom_to_modified;
        modified_base_atom_list(modified_nucleotide) = this_modified_base_atom_list;
        modified_base_to_hydrogens(modified_nucleotide) = this_modified_base_to_hydrogens;
        modified_base_to_hydrogen_coordinates(modified_nucleotide) = this_modified_base_to_hydrogen_coordinates;
    end


