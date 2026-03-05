import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

# =====================================================
# CONFIGURACIÓN DE PÁGINA
# =====================================================

st.set_page_config(
    page_title="Gaussian Input Generator",
    page_icon="⚛️"
)

st.sidebar.image("img/inGaussian.png", caption="Dr. Jesus Alvarado-Huayhuaz")

st.sidebar.title("¿Tienes las coordenadas?")
#st.sidebar.title("Gaussian Input Builder")


# =====================================================
# MOLÉCULA POR DEFECTO
# =====================================================

default_xyz = """3
Water molecule
O   0.0000   0.0000   0.0626
H  -0.7920   0.0000  -0.4973
H   0.7920   0.0000  -0.4973
"""

# =====================================================
# VISUALIZACIÓN 3D
# =====================================================

def show_xyz(xyz_data):
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(xyz_data, 'xyz')
    viewer.setStyle({'stick':{}})
    viewer.zoomTo()
    return viewer


# =====================================================
# SMILES → XYZ
# =====================================================

def smiles_to_xyz(smiles):

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    atoms = mol.GetAtoms()

    xyz_lines = [str(mol.GetNumAtoms()), "Generated from SMILES"]

    for atom in atoms:

        pos = conf.GetAtomPosition(atom.GetIdx())

        xyz_lines.append(
            f"{atom.GetSymbol():<2} {pos.x:>12.5f} {pos.y:>12.5f} {pos.z:>12.5f}"
        )

    return "\n".join(xyz_lines)


# =====================================================
# GENERADOR DE INPUT GAUSSIAN
# =====================================================

def convert_xyz_to_gaussian(
        xyz_content,
        method,
        basis,
        keywords,
        charge,
        multiplicity,
        memory,
        nproc,
        title,
        chk_name,
        print_orbitals,
        generate_cube):

    lines = xyz_content.strip().splitlines()
    atom_lines = lines[2:]

    keyword_list = keywords.copy()

    if print_orbitals:
        keyword_list.append("Pop=Full")

    if generate_cube:
        keyword_list.append("Density=Current")

    keyword_string = " ".join(keyword_list)

    route = f"# {method}/{basis} {keyword_string}"

    extra_section = ""
    
    if print_orbitals:
        extra_section += "\nPop=Full GFInput\n"
    
    if generate_cube:
        extra_section += "\nDensity=Current\n"
    
    header = f"""%chk={chk_name}

%mem={memory}
%nprocshared={nproc}
{route}

{title}

{charge} {multiplicity}
"""

    coord_block = "\n".join(atom_lines)

    cube_block = ""

    if generate_cube:
        cube_block = """

--Link1--
%chk={chk}
# {method}/{basis} Geom=AllCheck Guess=Read Density=Current

MEP Cube generation

{charge} {mult}

""".format(
            chk=chk_name,
            method=method,
            basis=basis,
            charge=charge,
            mult=multiplicity
        )

    return header + coord_block + "\n" + cube_block


# =====================================================
# PARÁMETROS
# =====================================================

def parameter_section():

    st.subheader("Parámetros de cálculo")

    col1, col2 = st.columns(2)

    with col1:

        method = st.text_input("Método", "B3LYP")
        basis = st.text_input("Conjunto Base", "6-31G(d)")

        keywords = st.multiselect(
            "Keywords de cálculo",
            ["Opt", "Freq", "TD", "NMR", "SCRF"],
            default=["Opt","Freq"]
        )

        charge = st.number_input("Carga", value=0, step=1)

    with col2:

        multiplicity = st.number_input("Multiplicidad", value=1, step=1)

        memory = st.text_input("Memoria", "2GB")

        nproc = st.number_input("Número de procesadores", value=8, step=1)

        chk_name = st.text_input("Nombre checkpoint (.chk)", "calculation.chk")

    title = st.text_input("Título del cálculo", "Gaussian Job")

    st.subheader("Opciones avanzadas")

    print_orbitals = st.checkbox("Imprimir orbitales moleculares")

    generate_cube = st.checkbox("Generar MEP cube para GaussView")

    return (
        method,
        basis,
        keywords,
        charge,
        multiplicity,
        memory,
        nproc,
        title,
        chk_name,
        print_orbitals,
        generate_cube
    )


# =====================================================
# MENÚ LATERAL
# =====================================================

menu = st.sidebar.radio(
    "¿Tienes coordenadas XYZ?",
    ("SI", "NO (usar SMILES)")
)

st.title("Generador de Input para Gaussian")

# =====================================================
# OPCIÓN 1 — XYZ
# =====================================================

if menu == "SI":

    uploaded_file = st.file_uploader("Subir archivo XYZ", type=["xyz"])

    if uploaded_file is not None:
        xyz_text = uploaded_file.read().decode("utf-8")
    else:
        st.info("Ejemplo: molécula de agua")
        xyz_text = default_xyz

    xyz_text = st.text_area("Contenido XYZ", xyz_text, height=200)

    st.subheader("Visualización 3D")

    viewer = show_xyz(xyz_text)

    st.components.v1.html(viewer._make_html(), height=400)

    (
        method,
        basis,
        keywords,
        charge,
        multiplicity,
        memory,
        nproc,
        title,
        chk_name,
        print_orbitals,
        generate_cube
    ) = parameter_section()

    gaussian_text = convert_xyz_to_gaussian(
        xyz_text,
        method,
        basis,
        keywords,
        charge,
        multiplicity,
        memory,
        nproc,
        title,
        chk_name,
        print_orbitals,
        generate_cube
    )

    st.subheader("Archivo Gaussian generado")

    st.text_area("Vista previa del archivo .com", gaussian_text, height=350)

    st.download_button(
        label="Descargar archivo .com",
        data=gaussian_text,
        file_name="gaussian_input.com",
        mime="text/plain"
    )


# =====================================================
# OPCIÓN 2 — SMILES
# =====================================================

else:

    smiles_input = st.text_input("Ingresa SMILES", "CCO")

    if smiles_input:

        try:

            xyz_text = smiles_to_xyz(smiles_input)

            st.subheader("Coordenadas generadas")

            st.text_area("XYZ generado", xyz_text, height=200)

            viewer = show_xyz(xyz_text)

            st.components.v1.html(viewer._make_html(), height=400)

            (
                method,
                basis,
                keywords,
                charge,
                multiplicity,
                memory,
                nproc,
                title,
                chk_name,
                print_orbitals,
                generate_cube
            ) = parameter_section()

            gaussian_text = convert_xyz_to_gaussian(
                xyz_text,
                method,
                basis,
                keywords,
                charge,
                multiplicity,
                memory,
                nproc,
                title,
                chk_name,
                print_orbitals,
                generate_cube
            )

            st.subheader("Archivo Gaussian generado")

            st.text_area("Vista previa del archivo .com", gaussian_text, height=350)

            st.download_button(
                label="Descargar archivo .com",
                data=gaussian_text,
                file_name="gaussian_input.com",
                mime="text/plain"
            )

        except:
            st.error("SMILES inválido")
