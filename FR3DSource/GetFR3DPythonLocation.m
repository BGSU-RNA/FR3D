% GetFR3DPythonLocation.m specifies where python code to read CIF files is located

if ~isempty(strfind(pwd,'zirbel')),
  PythonLocation = 'C:\Users\zirbel\Documents\GitHub\fr3d-python\examples\';
  PythonVersion = 'python';                    % python 2.7 required
else
  PythonLocation = '';
  PythonVersion = 'python';
end

