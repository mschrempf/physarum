use std::f64::consts::PI;

use pixels::{Pixels, SurfaceTexture};
use rand::Rng;
use rayon::prelude::*;
use winit::{
    dpi::PhysicalSize,
    event_loop::{ControlFlow, EventLoopBuilder},
    window::WindowBuilder,
};
use winit_input_helper::WinitInputHelper;

const GRID_SIZE: (usize, usize) = (500, 500);
const SENSE_DISTANCE: f64 = 50.0;
const SENSE_ANGLE: f64 = PI / 10.0;
const SPEED: f64 = 1.0;
const NOF_PARTICLES: usize = 1000;
const DEPOSIT: i32 = 10000;

struct Particle {
    x: f64,
    y: f64,
    heading: f64,
}

impl Particle {
    fn rotate(&mut self, grid: &IntensityGrid) {
        let mut rng = rand::thread_rng();
        let intensities = [-SENSE_ANGLE, 0.0, SENSE_ANGLE]
            .iter()
            .map(|angle_offset| {
                let x = self.x + (self.heading + angle_offset).cos() * SENSE_DISTANCE;
                let y = self.y + (self.heading + angle_offset).sin() * SENSE_DISTANCE;
                let mut sum = 0;
                for xd in -0..=0 {
                    for yd in -0..=0 {
                        sum += grid.at(x + xd as f64, y + yd as f64);
                    }
                }
                sum
            })
            .collect::<Vec<_>>();

        // no change
        if (intensities[1] > intensities[0]) && (intensities[1] > intensities[2]) {
            self.heading += rng.gen_range(-SENSE_ANGLE..=SENSE_ANGLE);
            return;
        }

        // random turn
        if (intensities[1] < intensities[0]) && (intensities[1] < intensities[2]) {
            self.heading += rng.gen_range(-SENSE_ANGLE..=SENSE_ANGLE);
            return;
        }

        // left turn
        if intensities[0] > intensities[2] {
            self.heading -= rng.gen_range(0.0..=0.5);
        } else {
            // right turn
            self.heading += rng.gen_range(0.0..=0.5);
        }
    }

    fn move_forward(&mut self) {
        self.x += SPEED * self.heading.cos();
        self.y += SPEED * self.heading.sin();

        // randomly bounce at the wall if needed
        let x_range = 0.0..(GRID_SIZE.0 as f64 - 1.0);
        let y_range = 0.0..(GRID_SIZE.1 as f64 - 1.0);

        if self.x < 0.0 {
            self.x = GRID_SIZE.0 as f64 - 1.0;
        } else if self.x > GRID_SIZE.0 as f64 - 1.0 {
            self.x = 0.0;
        }
        if self.y < 0.0 {
            self.y = GRID_SIZE.1 as f64 - 1.0;
        } else if self.y > GRID_SIZE.1 as f64 - 1.0 {
            self.y = 0.0;
        }

        //if !x_range.contains(&self.x) || !y_range.contains(&self.y) {
        //    self.x = self.x.clamp(x_range.start, x_range.end);
        //    self.y = self.y.clamp(y_range.start, y_range.end);
        //    self.heading = rand::thread_rng().gen_range(0.0..(2.0 * PI));
        //}
    }
}

struct IntensityGrid {
    data: [[i32; GRID_SIZE.0]; GRID_SIZE.1],
}

impl IntensityGrid {
    fn decay(&mut self) {
        self.data.par_iter_mut().for_each(|row| {
            row.iter_mut().for_each(|cell| {
                if *cell > 0 {
                    *cell -= 1;
                }
            })
        });
    }

    fn deposit(&mut self, particles: &[Particle]) {
        particles.iter().for_each(|particle| {
            let x = particle.x as usize;
            let y = particle.y as usize;

            self.data[x][y] = DEPOSIT;
        })
    }

    fn diffuse(&mut self) {
        let mut tmp = [[0i32; GRID_SIZE.0]; GRID_SIZE.1];

        for (x, row) in self.data.iter().enumerate() {
            for (y, _) in row.iter().enumerate() {
                let mut sum = 0;
                let mut nof_cells = 0;
                for xd in -1..=1 {
                    for yd in -1..=1 {
                        let xn = (x as i32) + xd;
                        let yn = (y as i32) + yd;

                        if (0..GRID_SIZE.0 as i32).contains(&xn)
                            && (0..GRID_SIZE.1 as i32).contains(&yn)
                        {
                            sum += self.data[xn as usize][yn as usize] as i32;
                            nof_cells += 1
                        }
                    }
                }

                tmp[x][y] = sum / nof_cells;
            }
        }

        self.data.copy_from_slice(&tmp);
    }

    fn at(&self, x: f64, y: f64) -> i32 {
        let x = x as usize;
        let y = y as usize;

        let x = x.clamp(0, GRID_SIZE.0 - 1);
        let y = y.clamp(0, GRID_SIZE.1 - 1);
        self.data[x][y]
    }
}

fn clear(frame: &mut [u8]) {
    for pixel in frame {
        *pixel = 0;
    }
}

fn draw_particles(particles: &[Particle], frame: &mut [u8]) {
    for particle in particles {
        let pixel_index: usize = 4 * (GRID_SIZE.0 * (particle.y as usize) + (particle.x as usize));

        let pixel = &mut frame[pixel_index..(pixel_index + 4)];
        pixel.copy_from_slice(&[255, 0, 0, 255]);
    }
}

fn draw_grid(grid: &IntensityGrid, frame: &mut [u8]) {
    for x in 0..GRID_SIZE.0 {
        for y in 0..GRID_SIZE.1 {
            let pixel_index: usize = 4 * (GRID_SIZE.0 * y + x);
            let pixel = &mut frame[pixel_index..(pixel_index + 4)];
            let tmp = ((grid.data[x][y] * u8::MAX as i32) / DEPOSIT).clamp(0, u8::MAX as i32) as u8;
            pixel.copy_from_slice(&[tmp; 4]);
        }
    }
}

fn main() {
    let mut grid = IntensityGrid {
        data: [[0; GRID_SIZE.0]; GRID_SIZE.1],
    };

    let mut rng = rand::thread_rng();
    let mut particles = (0..NOF_PARTICLES)
        .map(|_| Particle {
            x: rng.gen_range(0.0..((GRID_SIZE.0 - 1) as f64)),
            y: rng.gen_range(0.0..((GRID_SIZE.1 - 1) as f64)),
            heading: rng.gen_range(0.0..(2.0 * PI)),
        })
        .collect::<Vec<_>>();

    let event_loop = EventLoopBuilder::new().build();
    let window = WindowBuilder::new()
        .with_title("Slimey slimey :-)")
        .with_inner_size(PhysicalSize::new(GRID_SIZE.0 as u32, GRID_SIZE.1 as u32))
        .build(&event_loop)
        .unwrap();

    let mut input = WinitInputHelper::new();

    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
        Pixels::new(GRID_SIZE.0 as u32, GRID_SIZE.1 as u32, surface_texture).unwrap()
    };

    let mut diffuse_counter = 0;

    event_loop.run(move |event, _, control_flow| {
        // if let Event::RedrawRequested(_) = event {
        //     if let Err(_err) = pixels.render() {
        //         *control_flow = ControlFlow::ExitWithCode(-1);
        //         return;
        //     }
        // }

        // update world
        //println!("Update world ...");
        particles.par_iter_mut().for_each(|particle| {
            particle.rotate(&grid);
            particle.move_forward()
        });
        grid.deposit(&particles);
        grid.decay();
        if diffuse_counter == 0 {
            grid.diffuse();
        }
        diffuse_counter = (diffuse_counter + 1) % 2;
        //println!("... Finished");
        clear(pixels.frame_mut());
        draw_grid(&grid, pixels.frame_mut());
        pixels.render().unwrap();
        //draw_particles(&particles, pixels.frame_mut());

        //window.request_redraw();

        if input.update(&event) && input.close_requested() {
            *control_flow = ControlFlow::Exit;
        }
    });
}
