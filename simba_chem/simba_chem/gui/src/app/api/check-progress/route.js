export async function GET() {
  console.log('Progress check requested. Current progress:', global.solverProgress); // Debug log
  return Response.json({ 
    progress: global.solverProgress || 0 
  });
}